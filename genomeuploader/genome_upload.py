#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2017-2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import json
import logging
import os
import re
import sys
import xml.dom.minidom as minidom
import xml.etree.ElementTree as ET
from datetime import date
from datetime import datetime as dt

import ena
import pandas as pd
from constants import (
    BIN_MANDATORY_FIELDS,
    GEOGRAPHIC_LOCATIONS,
    GEOGRAPHY_DIGIT_COORDS,
    HQ,
    MAG_MANDATORY_FIELDS,
    METAGENOMES,
    MQ,
)
from ena import EnaQuery, EnaSubmit

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Create xmls and manifest files for genome upload and upload to ENA",
    )

    parser.add_argument("-u", "--upload_study", type=str, required=True, help="Study accession for genomes upload")
    parser.add_argument("--genome_info", type=str, required=True, help="Genomes metadata file")

    genomeType = parser.add_mutually_exclusive_group(required=True)
    genomeType.add_argument("-m", "--mags", action="store_true", help="Select for MAG upload")
    genomeType.add_argument("-b", "--bins", action="store_true", help="Select for bin upload")

    parser.add_argument("--out", type=str, help="Output folder. Default: working directory")
    parser.add_argument("--force", action="store_true", required=False, default=False, help="Forces reset of sample xml's backups")
    parser.add_argument(
        "--live",
        action="store_true",
        required=False,
        default=False,
        help="Uploads on ENA. Omitting this " + "option allows to validate samples beforehand",
    )
    parser.add_argument("--tpa", action="store_true", required=False, default=False, help="Select if uploading TPA-generated genomes")

    # Users can provide their credentials in environment variables or using a config file
    parser.add_argument("--centre_name", required=True, help="Name of the centre uploading genomes")
    parser.add_argument(
        "--private",
        required=False,
        action="store_true",
        default=False,
        help="if data is private",
    )

    return parser.parse_args(argv)


"""
Input table: expects the following parameters:
    genome_name: genome file name
    accessions: run(s) or assembly(ies) the genome was generated from
    assembly_software: assembler_vX.X
    binning_software: binner_vX.X
    binning_parameters: binning parameters
    stats_generation_software: software_vX.X
    completeness: float
    contamination: float
    rRNA_presence: True/False if 5S, 16S, and 23S genes have been detected in the genome
    NCBI_lineage: full NCBI lineage, either in tax id or strings
    broad_environment: string
    local_environment: string
    environmental_medium: string
    metagenome: string
    co-assembly: True/False, whether the genome was generated from a co-assembly
    genome_coverage : genome coverage
    genome_path: path to genome to upload
"""


def round_stats(stats):
    new_stat = round(float(stats), 2)
    if new_stat == 100.0:
        new_stat = 100
    if new_stat == 0:
        new_stat = 0.0

    return new_stat


def compute_mag_quality(completeness, contamination, rna_presence):
    rna_presence = str(rna_presence).lower() in ["true", "yes", "y"]
    quality = MQ
    if float(completeness) >= 90 and float(contamination) <= 5 and rna_presence:
        quality = HQ

    return quality, completeness, contamination


def extract_tax_info(tax_info):
    # if unclassified, block the execution
    lineage, kingdom_position_lineage, digit_annotation = tax_info.split(";"), 0, False
    lineage_first = lineage[0]
    if "Unclassified " in lineage_first:
        if "Archaea" in lineage_first:
            scientific_name = "uncultured archaeon"
        elif "Bacteria" in lineage_first:
            scientific_name = "uncultured bacterium"
        elif "Eukaryota" in lineage_first:
            scientific_name = "uncultured eukaryote"
        submittable, taxid, rank = ena.query_scientific_name(scientific_name, search_rank=True)
        return taxid, scientific_name

    kingdoms = ["Archaea", "Bacteria", "Eukaryota"]
    kingdom_taxa = ["2157", "2", "2759"]

    selected_kingdom, final_kingdom = kingdoms, ""
    if lineage[1].isdigit():
        selected_kingdom = kingdom_taxa
        kingdom_position_lineage = 2
        digit_annotation = True
    for index, k in enumerate(selected_kingdom):
        if digit_annotation:
            if k == lineage[kingdom_position_lineage]:
                final_kingdom = selected_kingdom[index]
                break
        else:
            if k in lineage[kingdom_position_lineage]:
                final_kingdom = selected_kingdom[index]
                break

    iterator = len(lineage) - 1
    submittable = False
    rank = ""
    while iterator != -1 and not submittable:
        scientific_name = lineage[iterator].strip()
        if digit_annotation:
            if "*" not in scientific_name:
                scientific_name = ena.query_taxid(scientific_name)
            else:
                iterator -= 1
                continue
        elif "__" in scientific_name:
            scientific_name = scientific_name.split("__")[1]
        else:
            raise ValueError("Unrecognised taxonomy format: " + scientific_name)
        submittable, taxid, rank = ena.query_scientific_name(scientific_name, search_rank=True)

        if not submittable:
            if final_kingdom == "Archaea" or final_kingdom == "2157":
                submittable, scientific_name, taxid = extract_archaea_info(scientific_name, rank)
            elif final_kingdom == "Bacteria" or final_kingdom == "2":
                submittable, scientific_name, taxid = extract_bacteria_info(scientific_name, rank)
            elif final_kingdom == "Eukaryota" or final_kingdom == "2759":
                submittable, scientific_name, taxid = extract_eukaryota_info(scientific_name, rank)
        iterator -= 1

    return taxid, scientific_name


def extract_eukaryota_info(name, rank):
    non_submittable = (False, "", 0)

    # Asterisks in given taxonomy suggest the classification might be not confident enough.
    if "*" in name:
        return non_submittable

    if rank == "super kingdom":
        name = "uncultured eukaryote"
        submittable, taxid = ena.query_scientific_name(name)
        return submittable, name, taxid
    else:
        name = name.capitalize() + " sp."
        submittable, taxid = ena.query_scientific_name(name)
        if submittable:
            return submittable, name, taxid
        else:
            name = "uncultured " + name
            submittable, taxid = ena.query_scientific_name(name)
            if submittable:
                return submittable, name, taxid
            else:
                name = name.replace(" sp.", "")
                submittable, taxid = ena.query_scientific_name(name)
                if submittable:
                    return submittable, name, taxid
                else:
                    return non_submittable


def extract_bacteria_info(name, rank):
    if rank == "species":
        name = name
    elif rank == "superkingdom":
        name = "uncultured bacterium".format()
    elif rank in ["family", "order", "class", "phylum"]:
        name = "uncultured {} bacterium".format(name)
    elif rank == "genus":
        name = "uncultured {} sp.".format(name)

    submittable, taxid, rank = ena.query_scientific_name(name, search_rank=True)
    if not submittable:
        if rank in ["species", "genus"] and name.lower().endswith("bacteria"):
            name = "uncultured {}".format(name.lower().replace("bacteria", "bacterium"))
        elif rank == "family":
            if name.lower() == "deltaproteobacteria":
                name = "uncultured delta proteobacterium"
        submittable, taxid = ena.query_scientific_name(name)

    return submittable, name, taxid


def extract_archaea_info(name, rank):
    if rank == "species":
        name = name
    elif rank == "superkingdom":
        name = "uncultured archaeon"
    elif rank == "phylum":
        if "Euryarchaeota" in name:
            name = "uncultured euryarchaeote"
        elif "Candidatus" in name:
            name = "{} archaeon".format(name)
        else:
            name = "uncultured {} archaeon".format(name)
    elif rank in ["family", "order", "class"]:
        name = "uncultured {} archaeon".format(name)
    elif rank == "genus":
        name = "uncultured {} sp.".format(name)

    submittable, taxid, rank = ena.query_scientific_name(name, search_rank=True)
    if not submittable:
        if "Candidatus" in name:
            if rank == "phylum":
                name = name.replace("Candidatus ", "")
            elif rank == "family":
                name = name.replace("uncultured ", "")
            submittable, taxid = ena.query_scientific_name(name)

    return submittable, name, taxid


def multiple_element_set(metadata_list):
    return len(set(metadata_list)) > 1


def combine_ena_info(genome_info, ena_dict):
    for g in genome_info:
        # TODO: optimise all the part below
        if genome_info[g]["co-assembly"]:
            instrument_list, collection_list, country_list = [], [], []
            study_list, description_list, samples_list = [], [], []
            long_list, latit_list = [], []
            for run in genome_info[g]["accessions"]:
                instrument_list.append(ena_dict[run]["instrumentModel"])
                collection_list.append(ena_dict[run]["collectionDate"])
                country_list.append(ena_dict[run]["country"])
                study_list.append(ena_dict[run]["study"])
                description_list.append(ena_dict[run]["projectDescription"])
                samples_list.append(ena_dict[run]["sampleAccession"])
                long_list.append(ena_dict[run]["longitude"])
                latit_list.append(ena_dict[run]["latitude"])

            genome_info[g]["study"] = study_list[0]
            genome_info[g]["description"] = description_list[0]

            instrument = instrument_list[0]
            if multiple_element_set(instrument_list):
                instrument = ",".join(instrument_list)
            genome_info[g]["sequencingMethod"] = instrument

            collection_date = collection_list[0]
            if multiple_element_set(collection_list):
                collection_date = "not provided"
            genome_info[g]["collectionDate"] = collection_date

            country = country_list[0]
            if multiple_element_set(country_list):
                country = "not applicable"
            genome_info[g]["country"] = country

            latitude = latit_list[0]
            if multiple_element_set(latit_list):
                latitude = "not provided"
            genome_info[g]["latitude"] = str(round(float(latitude), GEOGRAPHY_DIGIT_COORDS))

            longitude = long_list[0]
            if multiple_element_set(long_list):
                longitude = "not provided"
            genome_info[g]["longitude"] = str(round(float(longitude), GEOGRAPHY_DIGIT_COORDS))

            samples = samples_list[0]
            if multiple_element_set(samples_list):
                samples = ",".join(samples_list)
            genome_info[g]["sample_accessions"] = samples
        else:
            run = genome_info[g]["accessions"][0]
            genome_info[g]["sequencingMethod"] = ena_dict[run]["instrumentModel"]
            genome_info[g]["collectionDate"] = ena_dict[run]["collectionDate"]
            genome_info[g]["study"] = ena_dict[run]["study"]
            genome_info[g]["description"] = ena_dict[run]["projectDescription"]
            genome_info[g]["sample_accessions"] = ena_dict[run]["sampleAccession"]
            genome_info[g]["country"] = ena_dict[run]["country"]
            genome_info[g]["longitude"] = ena_dict[run]["longitude"]
            genome_info[g]["latitude"] = ena_dict[run]["latitude"]

        genome_info[g]["accessions"] = ",".join(genome_info[g]["accessions"])


def get_accessions(accessions_file):
    accession_dict = {}
    with open(accessions_file, "r") as f:
        for line in f:
            line = line.split("\t")
            alias = line[0]
            accession = line[1].rstrip("\n")
            accession_dict[alias] = accession

    return accession_dict


def save_accessions(alias_accession_dict, accessions_file, write_mode):
    with open(accessions_file, write_mode) as f:
        for elem in alias_accession_dict:
            f.write("{}\t{}\n".format(elem, alias_accession_dict[elem]))


def create_manifest_dictionary(run, alias, assembly_software, sequencing_method, mag_path, gen, study, coverage, is_coassembly):
    manifest_dict = {
        "accessions": run,
        "alias": alias,
        "assembler": assembly_software,
        "sequencingMethod": sequencing_method,
        "genome_path": mag_path,
        "genome_name": gen,
        "study": study,
        "coverageDepth": coverage,
        "co-assembly": is_coassembly,
    }

    return manifest_dict


def compute_manifests(genomes):
    manifest_info = {}
    for g in genomes:
        manifest_info[g] = create_manifest_dictionary(
            genomes[g]["accessions"],
            genomes[g]["alias"],
            genomes[g]["assembly_software"],
            genomes[g]["sequencingMethod"],
            genomes[g]["genome_path"],
            g,
            genomes[g]["study"],
            genomes[g]["genome_coverage"],
            genomes[g]["co-assembly"],
        )

    return manifest_info


def get_study_from_xml(sample):
    description = sample.childNodes[5].childNodes[0].data
    study = description.split(" ")[-1][:-1]

    return study


def recover_info_from_xml(genome_dict, sample_xml, live_mode):
    logger.info("Retrieving data for genome submission...")

    # extract list of genomes (samples) to be registered
    xml_structure = minidom.parse(sample_xml)
    samples = xml_structure.getElementsByTagName("SAMPLE")

    for s in samples:
        study = get_study_from_xml(s)

        # extract alias from xml and find a match with genomes the user is uploading
        xml_alias = s.attributes["alias"].value
        if not live_mode:  # remove time stamp if test mode is selected
            alias_split = xml_alias.split("_")
            xml_alias = "_".join(alias_split[:-1])
        for gen in genome_dict:
            # if match is found, associate attributes listed in the xml file
            # with genomes to upload
            if xml_alias == gen:
                if not live_mode:
                    current_time_stamp = str(int(dt.timestamp(dt.now())))
                    xml_alias = gen + "_" + current_time_stamp
                    s.attributes["alias"].value = xml_alias
                    sample_title = s.getElementsByTagName("TITLE")[0]
                    sample_title_value = sample_title.firstChild.nodeValue.split("_")
                    sample_title_value[-1] = current_time_stamp
                    new_sample_title = "_".join(sample_title_value)
                    s.getElementsByTagName("TITLE")[0].firstChild.replaceWholeText(new_sample_title)
                attributes = s.childNodes[7].getElementsByTagName("SAMPLE_ATTRIBUTE")
                seq_method, ass_software = "", ""
                for a in attributes:
                    tag_elem = a.getElementsByTagName("TAG")
                    tag = tag_elem[0].childNodes[0].nodeValue
                    if tag == "sequencing method":
                        seq_method_elem = a.getElementsByTagName("VALUE")
                        seq_method = seq_method_elem[0].childNodes[0].nodeValue
                    elif tag == "assembly software":
                        ass_software_elem = a.getElementsByTagName("VALUE")
                        ass_software = ass_software_elem[0].childNodes[0].nodeValue
                    if not seq_method == "" and not ass_software == "":
                        break

                genome_dict[gen]["accessions"] = ",".join(genome_dict[gen]["accessions"])
                genome_dict[gen]["alias"] = xml_alias
                genome_dict[gen]["assembly_software"] = ass_software
                genome_dict[gen]["sequencingMethod"] = seq_method
                genome_dict[gen]["study"] = study
                break

    if not live_mode:
        for s in samples:
            xml_structure.firstChild.appendChild(s)

        with open(sample_xml, "wb") as f:
            dom_string = xml_structure.toprettyxml().encode("utf-8")
            dom_string = b"\n".join([s for s in dom_string.splitlines() if s.strip()])
            f.write(dom_string)


def create_sample_attribute(sample_attributes, data_list, mag_data=None):
    tag = data_list[0]
    value = data_list[1]
    if mag_data:
        value = str(mag_data[value])
    units = None
    if len(data_list) == 3:
        units = data_list[2]

    new_sample_attr = ET.SubElement(sample_attributes, "SAMPLE_ATTRIBUTE")
    ET.SubElement(new_sample_attr, "TAG").text = tag
    ET.SubElement(new_sample_attr, "VALUE").text = value

    if units:
        ET.SubElement(new_sample_attr, "UNITS").text = units


def write_genomes_xml(genomes, xml_path, genome_type, centre_name, tpa):
    map_sample_attributes = [
        # tag - value - unit (optional)
        ["project name", "description"],
        ["sequencing method", "sequencingMethod"],
        ["assembly software", "assembly_software"],
        ["assembly quality", "MAG_quality"],
        ["binning software", "binning_software"],
        ["binning parameters", "binning_parameters"],
        ["completeness software", "stats_generation_software"],
        ["completeness score", "completeness", "%"],
        ["contamination score", "contamination", "%"],
        ["isolation source", "isolationSource"],
        ["collection date", "collectionDate"],
        ["geographic location (country and/or sea)", "country"],
        ["geographic location (latitude)", "latitude", "DD"],
        ["geographic location (longitude)", "longitude", "DD"],
        ["broad-scale environmental context", "broad_environment"],
        ["local environmental context", "local_environment"],
        ["environmental medium", "environmental_medium"],
        ["sample derived from", "sample_accessions"],
        ["metagenomic source", "metagenome"],
    ]

    checklist, assembly_type = "ERC000047", "Metagenome-assembled genome"
    if genome_type == "bins":
        checklist = "ERC000050"
        assembly_type = "binned metagenome"

    constant_sample_attributes = [
        # tag - value
        ["taxonomic identity marker", "multi-marker approach"],
        ["investigation type", "metagenome-assembled genome"],
        ["ENA-CHECKLIST", checklist],
    ]

    tpa_description = ""
    if tpa:
        tpa_description = "Third Party Annotation (TPA) "

    sample_set = ET.Element("SAMPLE_SET")

    for g in genomes:
        plural = ""
        if genomes[g]["co-assembly"]:
            plural = "s"
        description = "This sample represents a {}{} assembled from the " "metagenomic run{} {} of study {}.".format(
            tpa_description, assembly_type, plural, genomes[g]["accessions"], genomes[g]["study"]
        )

        sample = ET.SubElement(sample_set, "SAMPLE")
        sample.set("alias", genomes[g]["alias"])
        sample.set("center_name", centre_name)

        ET.SubElement(sample, "TITLE").text = "{}: {}".format(assembly_type, genomes[g]["alias"])
        sample_name = ET.SubElement(sample, "SAMPLE_NAME")
        ET.SubElement(sample_name, "TAXON_ID").text = genomes[g]["taxID"]
        ET.SubElement(sample_name, "SCIENTIFIC_NAME").text = genomes[g]["scientific_name"]
        ET.SubElement(sample_name, "COMMON_NAME")

        ET.SubElement(sample, "DESCRIPTION").text = description

        sample_attributes = ET.SubElement(sample, "SAMPLE_ATTRIBUTES")

        for mapping in map_sample_attributes:
            create_sample_attribute(sample_attributes, mapping, genomes[g])

        for constant in constant_sample_attributes:
            create_sample_attribute(sample_attributes, constant)

    with open(xml_path, "wb") as f:
        dom = minidom.parseString(ET.tostring(sample_set, encoding="utf-8"))
        f.write(dom.toprettyxml().encode("utf-8"))


def write_submission_xml(upload_dir, centre_name, study=True):
    today = str(date.today())
    sub_xml = os.path.join(upload_dir, "submission.xml")

    submission = ET.Element("SUBMISSION")
    submission.set("center_name", centre_name)

    # template
    actions = ET.SubElement(submission, "ACTIONS")
    action_sub = ET.SubElement(actions, "ACTION")
    ET.SubElement(action_sub, "ADD")

    # attributes: function and hold date
    if study:
        action_hold = ET.SubElement(actions, "ACTION")
        hold = ET.SubElement(action_hold, "HOLD")
        hold.set("HoldUntilDate", today)

    with open(sub_xml, "wb") as submission_file:
        dom = minidom.parseString(ET.tostring(submission, encoding="utf-8"))
        submission_file.write(dom.toprettyxml().encode("utf-8"))

    return sub_xml


def generate_genome_manifest(genome_info, study, manifests_root, alias_to_sample, genome_type, tpa):
    manifest_path = os.path.join(manifests_root, f'{genome_info["genome_name"]}.manifest')

    tpa_addition, multiple_runs = "", ""
    if tpa:
        tpa_addition = "Third Party Annotation (TPA) "
    if genome_info["co-assembly"]:
        multiple_runs = "s"
    assembly_type = "Metagenome-Assembled Genome (MAG)"
    if genome_type == "bins":
        assembly_type = "binned metagenome"

    values = (
        ("STUDY", study),
        ("SAMPLE", alias_to_sample[genome_info["alias"]]),
        ("ASSEMBLYNAME", genome_info["alias"]),
        ("ASSEMBLY_TYPE", assembly_type),
        ("COVERAGE", genome_info["coverageDepth"]),
        ("PROGRAM", genome_info["assembler"]),
        ("PLATFORM", genome_info["sequencingMethod"]),
        ("MOLECULETYPE", "genomic DNA"),
        (
            "DESCRIPTION",
            (
                "This is a {}bin derived from the primary whole genome "
                "shotgun (WGS) data set {}. This sample represents a {} from the "
                "metagenomic run{} {}.".format(tpa_addition, genome_info["study"], assembly_type, multiple_runs, genome_info["accessions"])
            ),
        ),
        ("RUN_REF", genome_info["accessions"]),
        ("FASTA", os.path.abspath(genome_info["genome_path"])),
    )
    logger.info("Writing manifest file (.manifest) for {}.".format(genome_info["alias"]))
    with open(manifest_path, "w") as outfile:
        for k, v in values:
            manifest = f"{k}\t{v}\n"
            outfile.write(manifest)
        if tpa:
            outfile.write("TPA\ttrue\n")


class GenomeUpload:
    def __init__(
        self,
        upload_study: str,
        centre_name: str,
        genome_info: str,
        bins: bool = False,
        live: bool = False,
        private: bool = False,
        tpa: bool = False,
        force: bool = False,
        out: str = None,
    ):
        """
        Submission of genomes.

        :param upload_study: Study accession for genomes upload.
        :param centre_name: Name of the centre uploading genomes.
        :param genome_info: Genomes metadata file.
        :param bins: Performs bin upload.
        :params live: Live upload to ENA.
        :param private: Is this a private study?
        :param tpa: Is this a third-party assembly?
        :param force: Resets sample XML backups.
        :param out: Output folder.
        """

        self.genome_type = "bins" if bins else "MAGs"
        self.live = live
        self.private = private

        self.tpa = tpa
        self.centre_name = centre_name
        self.force = force

        self.work_dir = out if out is not None else os.getcwd()
        self.upload_dir = self.generate_genomes_upload_dir()

        if not upload_study:
            raise ValueError("No project selected for genome upload [-u, --upload_study].")
        self.upload_study = upload_study

        if not os.path.exists(genome_info):
            raise FileNotFoundError(f"Genome metadata file {genome_info} does not exist")
        self.genome_metadata = genome_info

    def validate_metadata_tsv(self):
        logger.info("Retrieving info for genomes to submit...")

        all_fields = MAG_MANDATORY_FIELDS + BIN_MANDATORY_FIELDS
        metadata = pd.read_csv(self.genome_metadata, sep="\t", usecols=all_fields)

        # make sure there are no empty cells
        clean_columns = list(metadata.dropna(axis=1))
        if self.genome_type == "MAGs":
            missing_values = [item for item in all_fields if item not in clean_columns]
        else:
            missing_values = [item for item in BIN_MANDATORY_FIELDS if item not in clean_columns]

        if missing_values:
            raise ValueError(
                "The following mandatory fields have missing values in " + "the input file: {}".format(", ".join(missing_values))
            )

        # check amount of genomes to register at the same time
        if len(metadata) >= 5000:
            raise ValueError("Genomes need to be registered in batches of 5000 genomes or smaller.")

        # check whether accessions follow the right format
        accessions_reg_exp = re.compile(r"([E|S|D]R[R|Z]\d{6,})")

        accession_comparison = pd.DataFrame(columns=["genome_name", "attemptive_accessions", "correct", "mismatching", "co-assembly"])
        accession_comparison["genome_name"] = metadata["genome_name"]

        accession_comparison["attemptive_accessions"] = metadata["accessions"].map(lambda a: len(a.split(",")))

        accession_comparison["correct"] = metadata["accessions"].map(lambda a: len(accessions_reg_exp.findall(a)))

        accession_comparison["mismatching"] = accession_comparison.apply(
            lambda row: True if row["attemptive_accessions"] == row["correct"] else None, axis=1
        ).isna()

        mismatching_accessions = accession_comparison[accession_comparison["mismatching"]]["genome_name"]
        if not mismatching_accessions.empty:
            raise ValueError(
                "Run accessions are not correctly formatted for the following " + "genomes: " + ",".join(mismatching_accessions.values)
            )

        # check whether completeness and contamination are floats
        try:
            pd.to_numeric(metadata["completeness"])
            pd.to_numeric(metadata["contamination"])
            pd.to_numeric(metadata["genome_coverage"])
        except ValueError:
            raise ValueError("Completeness, contamination or coverage values should be formatted as floats")

        # check whether all co-assemblies have more than one run associated and viceversa
        accession_comparison["co-assembly"] = metadata["co-assembly"]
        coassembly_discrepancy = metadata[
            ((accession_comparison["correct"] < 2) & (accession_comparison["co-assembly"]))
            | ((accession_comparison["correct"] > 1) & (~accession_comparison["co-assembly"]))
        ]["genome_name"]
        if not coassembly_discrepancy.empty:
            raise ValueError(
                "The following genomes show discrepancy between number of runs "
                "involved and co-assembly status: " + ",".join(coassembly_discrepancy.values)
            )

        # are provided metagenomes part of the accepted metagenome list?
        if False in metadata.apply(lambda row: True if row["metagenome"] in METAGENOMES else False, axis=1).unique():
            raise ValueError("Metagenomes associated with each genome need to belong to ENA's " + "approved metagenomes list.")

        # do provided file paths exist?
        if False in metadata.apply(lambda row: True if os.path.exists(row["genome_path"]) else False, axis=1).unique():
            raise FileNotFoundError("Some genome paths do not exist.")

        # check genome name lengths
        # if not (metadata["genome_name"].map(lambda a: len(a) < 20).all()):
        #    raise ValueError("Genome names must be shorter than 20 characters.")

        # create dictionary while checking genome name uniqueness
        uniqueness = metadata["genome_name"].nunique() == metadata["genome_name"].size
        if uniqueness:
            if not self.live:
                timestamp = str(int(dt.timestamp(dt.now())))
                timestamp_names = [row["genome_name"] + "_" + timestamp for index, row in metadata.iterrows()]
                metadata["unique_genome_name"] = timestamp_names
                genome_info = metadata.set_index("unique_genome_name").transpose().to_dict()
            else:
                genome_info = metadata.set_index("genome_name").transpose().to_dict()
        else:
            raise ValueError("Duplicate names found in genome names")

        return genome_info

    def extract_genomes_info(self):
        genome_info = self.validate_metadata_tsv()
        for gen in genome_info:
            genome_info[gen]["accessions"] = genome_info[gen]["accessions"].split(",")
            accession_type = "run"
            assembly_reg_exp = re.compile(r"([E|S|D]RZ\d{6,})")
            if assembly_reg_exp.findall(genome_info[gen]["accessions"][0]):
                accession_type = "assembly"
            genome_info[gen]["accessionType"] = accession_type

            genome_info[gen]["isolationSource"] = genome_info[gen]["metagenome"]

            try:
                (
                    genome_info[gen]["MAG_quality"],
                    genome_info[gen]["completeness"],
                    genome_info[gen]["contamination"],
                ) = compute_mag_quality(
                    str(round_stats(genome_info[gen]["completeness"])),
                    str(round_stats(genome_info[gen]["contamination"])),
                    genome_info[gen]["rRNA_presence"],
                )
            except IndexError:
                pass

            if str(genome_info[gen]["co-assembly"]).lower() in ["yes", "y", "true"]:
                genome_info[gen]["co-assembly"] = True
            else:
                genome_info[gen]["co-assembly"] = False

            genome_info[gen]["alias"] = gen

            tax_id, scientific_name = extract_tax_info(genome_info[gen]["NCBI_lineage"])
            genome_info[gen]["taxID"] = tax_id
            genome_info[gen]["scientific_name"] = scientific_name

        return genome_info

    def extract_ena_info(self, genome_info):
        logger.info("Retrieving project and run info from ENA (this might take a while)...")

        # retrieving metadata from runs (and runs from assembly accessions if provided)
        all_runs = []
        for g in genome_info:
            if genome_info[g]["accessionType"] == "assembly":
                derived_runs = []
                for acc in genome_info[g]["accessions"]:
                    ena_query = EnaQuery(acc, "run_assembly", self.private)
                    derived_runs.append(ena_query.build_query())
                genome_info[g]["accessions"] = derived_runs
            all_runs.extend(genome_info[g]["accessions"])

        runs_set, study_set, samples_dict, temp_dict = set(all_runs), set(), {}, {}
        for r in runs_set:
            ena_query = EnaQuery(r, "run", self.private)
            run_info = ena_query.build_query()
            study_set.add(run_info["secondary_study_accession"])
            samples_dict[r] = run_info["sample_accession"]

        if not study_set:
            raise ValueError("No study corresponding to runs found.")

        backup_file = os.path.join(self.upload_dir, "ENA_backup.json")
        counter = 0
        if not os.path.exists(backup_file):
            with open(backup_file, "w") as file:
                pass
        with open(backup_file, "r+") as file:
            try:
                backup_dict = json.load(file)
                temp_dict = dict(backup_dict)
                logger.info(f"A backup file {backup_file} for ENA sample metadata has been found.")
            except json.decoder.JSONDecodeError:
                backup_dict = {}
            for s in study_set:
                ena_query = EnaQuery(s, "study", self.private)
                study_info = ena_query.build_query()
                project_description = study_info["study_description"]

                ena_query = EnaQuery(s, "study_runs", self.private)
                ena_info = ena_query.build_query()
                if ena_info == []:
                    raise IOError("No runs found on ENA for project {}.".format(s))

                for run, item in enumerate(ena_info):
                    run_accession = ena_info[run]["run_accession"]
                    if run_accession not in backup_dict:
                        if run_accession in runs_set:
                            sample_accession = ena_info[run]["sample_accession"]
                            ena_query = EnaQuery(sample_accession, "sample", self.private)
                            sample_info = ena_query.build_query()

                            if "latitude" in sample_info and "longitude" in sample_info:
                                latitude = sample_info["latitude"]
                                longitude = sample_info["longitude"]
                            else:
                                location = sample_info["location"]
                                latitude, longitude = None, None
                                if "N" in location:
                                    latitude = location.split("N")[0].strip()
                                    longitude = location.split("N")[1].strip()
                                elif "S" in location:
                                    latitude = "-" + location.split("S")[0].strip()
                                    longitude = location.split("S")[1].strip()

                                if "W" in longitude:
                                    longitude = "-" + longitude.split("W")[0].strip()
                                elif longitude.endswith("E"):
                                    longitude = longitude.split("E")[0].strip()

                            if latitude:
                                latitude = "{:.{}f}".format(round(float(latitude), GEOGRAPHY_DIGIT_COORDS), GEOGRAPHY_DIGIT_COORDS)
                            else:
                                latitude = "not provided"

                            if longitude:
                                longitude = "{:.{}f}".format(round(float(longitude), GEOGRAPHY_DIGIT_COORDS), GEOGRAPHY_DIGIT_COORDS)
                            else:
                                longitude = "not provided"

                            country = sample_info["country"].split(":")[0]
                            if country not in GEOGRAPHIC_LOCATIONS:
                                country = "not provided"

                            collectionDate = sample_info["collection_date"]
                            if collectionDate == "" or collectionDate == "missing":
                                collectionDate = "not provided"

                            temp_dict[run_accession] = {
                                "instrumentModel": ena_info[run]["instrument_model"],
                                "collectionDate": collectionDate,
                                "country": country,
                                "latitude": latitude,
                                "longitude": longitude,
                                "projectDescription": project_description,
                                "study": s,
                                "sampleAccession": samples_dict[run_accession],
                            }
                            counter += 1

                            if (counter % 10 == 0) or (len(runs_set) - len(backup_dict) == counter):
                                file.seek(0)
                                file.write(json.dumps(temp_dict))
                                file.truncate()
        temp_dict = {**temp_dict, **backup_dict}
        combine_ena_info(genome_info, temp_dict)

    def generate_genomes_upload_dir(self):
        upload_name = "MAG_upload"
        if self.genome_type == "bins":
            upload_name = upload_name.replace("MAG", "bin")
        upload_dir = os.path.join(self.work_dir, upload_name)
        os.makedirs(upload_dir, exist_ok=True)
        return upload_dir

    def create_genome_dictionary(self, samples_xml):
        logger.info("Retrieving data for MAG submission...")

        genome_info = self.extract_genomes_info()

        if not os.path.exists(samples_xml) or self.force:
            self.extract_ena_info(genome_info)
            logger.info("Writing genome registration XML...")

            write_genomes_xml(genome_info, samples_xml, self.genome_type, self.centre_name, self.tpa)
            logger.info("All files have been written to " + self.upload_dir)
        else:
            recover_info_from_xml(genome_info, samples_xml, self.live)

        return genome_info

    def genome_upload(self):
        if not self.live:
            logger.warning("Warning: genome submission is not in live mode, " + "files will be validated, but not uploaded.")
        samples_xml = os.path.join(self.upload_dir, "genome_samples.xml")
        submission_xml = os.path.join(self.upload_dir, "submission.xml")
        genomes, manifest_info = {}, {}

        # submission xml existence
        if not os.path.exists(submission_xml):
            submission_xml = write_submission_xml(self.upload_dir, self.centre_name, False)

        # sample xml generation or recovery
        genomes = self.create_genome_dictionary(samples_xml)

        # manifests creation
        manifest_dir = os.path.join(self.upload_dir, "manifests")
        os.makedirs(manifest_dir, exist_ok=True)

        accessionsgen = "registered_MAGs.tsv"
        if self.genome_type == "bins":
            accessionsgen = accessionsgen.replace("MAG", "bin")
        if not self.live:
            accessionsgen = accessionsgen.replace(".tsv", "_test.tsv")

        accessions_file = os.path.join(self.upload_dir, accessionsgen)
        save = False
        write_mode = "a"
        if os.path.exists(accessions_file):
            if not self.live:
                save = True
                if self.force:
                    write_mode = "w"
            if not save:
                logger.info("Genome samples already registered, reading ERS accessions...")
                alias_to_new_sample_accession = get_accessions(accessions_file)
        else:
            save = True

        if save:
            logger.info("Registering genome samples XMLs...")
            ena_submit = EnaSubmit(samples_xml, submission_xml, self.live)
            alias_to_new_sample_accession = ena_submit.handle_genomes_registration()
            save_accessions(alias_to_new_sample_accession, accessions_file, write_mode)

        logger.info("Generating manifest files...")

        manifest_info = compute_manifests(genomes)

        for m in manifest_info:
            generate_genome_manifest(
                manifest_info[m], self.upload_study, manifest_dir, alias_to_new_sample_accession, self.genome_type, self.tpa
            )


def main():
    args = parse_args(sys.argv[1:])
    ena_upload = GenomeUpload(
        upload_study=args.upload_study,
        centre_name=args.centre_name,
        genome_info=args.genome_info,
        bins=args.bins,
        live=args.live,
        private=args.private,
        tpa=args.tpa,
        force=args.force,
        out=args.out,
    )
    ena_upload.genome_upload()


if __name__ == "__main__":
    main()
    logger.info("Completed")
