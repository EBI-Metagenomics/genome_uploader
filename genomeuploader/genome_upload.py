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

import json
import logging
import importlib.metadata
import re
import click
import pandas as pd
from datetime import date, datetime as dt
from pathlib import Path

import xml.etree.ElementTree as et
import xml.dom.minidom as minidom

import genomeuploader.ena as ena
from genomeuploader.constants import (
    BIN_MANDATORY_FIELDS,
    GEOGRAPHIC_LOCATIONS,
    GEOGRAPHY_DIGIT_COORDS,
    HQ,
    MAG_MANDATORY_FIELDS,
    METAGENOMES,
    MQ,
)
from genomeuploader.ena import EnaQuery
from genomeuploader.ena_submit import EnaSubmit

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

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
    elif rank == "domain":
        name = "uncultured bacterium"
    elif rank in ["family", "order", "class", "phylum"]:
        name = f"uncultured {name} bacterium"
    elif rank == "genus":
        name = f"uncultured {name} sp."

    submittable, taxid, rank = ena.query_scientific_name(name, search_rank=True)
    if not submittable:
        if rank in ["species", "genus"] and name.lower().endswith("bacteria"):
            name = f"uncultured {name.lower().replace('bacteria', 'bacterium')}"
        elif rank == "family":
            if name.lower() == "deltaproteobacteria":
                name = "uncultured delta proteobacterium"
        submittable, taxid = ena.query_scientific_name(name)

    return submittable, name, taxid


def extract_archaea_info(name, rank):
    if rank == "species":
        name = name
    elif rank == "domain":
        name = "uncultured archaeon"
    elif rank == "phylum":
        if "Euryarchaeota" in name:
            name = "uncultured euryarchaeote"
        elif "Candidatus" in name:
            name = f"{name} archaeon"
        else:
            name = f"uncultured {name} archaeon"
    elif rank in ["family", "order", "class"]:
        name = f"uncultured {name} archaeon"
    elif rank == "genus":
        name = f"uncultured {name} sp."

    submittable, taxid, rank = ena.query_scientific_name(name, search_rank=True)
    if not submittable:
        if "Candidatus" in name:
            if rank == "phylum":
                name = name.replace("Candidatus ", "")
            elif rank == "family":
                name = name.replace("uncultured ", "")
            submittable, taxid = ena.query_scientific_name(name)

    return submittable, name, taxid

# NOT PRESENT IN MAIN
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
            if collection_date.lower() in ["not available", "na"]:
                collection_date = "missing: third party data"
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
            if ena_dict[run]["collectionDate"].lower() in ["not applicable", "not available", "na"]:
                genome_info[g]["collectionDate"] = "not provided"
            else:
                genome_info[g]["collectionDate"] = ena_dict[run]["collectionDate"]
            genome_info[g]["study"] = ena_dict[run]["study"]
            genome_info[g]["description"] = ena_dict[run]["projectDescription"]
            genome_info[g]["sample_accessions"] = ena_dict[run]["sampleAccession"]
            genome_info[g]["country"] = ena_dict[run]["country"]
            genome_info[g]["longitude"] = ena_dict[run]["longitude"]
            genome_info[g]["latitude"] = ena_dict[run]["latitude"]

        genome_info[g]["accessions"] = ",".join(genome_info[g]["accessions"])

def save_accessions(alias_accession_dict, accessions_file, write_mode):
    with accessions_file.open(write_mode) as f:
        for elem in alias_accession_dict:
            f.write(f"{elem}\t{alias_accession_dict[elem]}\n")


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

def create_sample_attribute(sample_attributes, data_list, mag_data=None):
    tag = data_list[0]
    value = data_list[1]
    if mag_data:
        value = str(mag_data[value])
    units = None
    if len(data_list) == 3:
        units = data_list[2]

    new_sample_attr = et.SubElement(sample_attributes, "SAMPLE_ATTRIBUTE")
    et.SubElement(new_sample_attr, "TAG").text = tag
    et.SubElement(new_sample_attr, "VALUE").text = value

    if units:
        et.SubElement(new_sample_attr, "UNITS").text = units


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
        ["isolation_source", "isolationSource"],
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

    sample_set = et.Element("SAMPLE_SET")

    for g in genomes:
        plural = ""
        if genomes[g]["co-assembly"]:
            plural = "s"
        description = f"This sample represents a {tpa_description}{assembly_type} assembled from the metagenomic run{plural} {genomes[g]['accessions']} of study {genomes[g]['study']}."

        sample = et.SubElement(sample_set, "SAMPLE")
        sample.set("alias", genomes[g]["alias"])
        sample.set("center_name", centre_name)

        et.SubElement(sample, "TITLE").text = f"{assembly_type}: {genomes[g]['alias']}"
        sample_name = et.SubElement(sample, "SAMPLE_NAME")
        et.SubElement(sample_name, "TAXON_ID").text = genomes[g]["taxID"]
        et.SubElement(sample_name, "SCIENTIFIC_NAME").text = genomes[g]["scientific_name"]
        et.SubElement(sample_name, "COMMON_NAME")

        et.SubElement(sample, "DESCRIPTION").text = description

        sample_attributes = et.SubElement(sample, "SAMPLE_ATTRIBUTES")

        for mapping in map_sample_attributes:
            create_sample_attribute(sample_attributes, mapping, genomes[g])

        for constant in constant_sample_attributes:
            create_sample_attribute(sample_attributes, constant)

    with xml_path.open("wb") as f:
        dom = minidom.parseString(et.tostring(sample_set, encoding="utf-8"))
        f.write(dom.toprettyxml().encode("utf-8"))


def write_submission_xml(upload_dir, centre_name, study=True):
    today = str(date.today())
    sub_xml = upload_dir / "submission.xml"

    submission = et.Element("SUBMISSION")
    submission.set("center_name", centre_name)

    # template
    actions = et.SubElement(submission, "ACTIONS")
    action_sub = et.SubElement(actions, "ACTION")
    et.SubElement(action_sub, "ADD")

    # attributes: function and hold date
    if study:
        action_hold = et.SubElement(actions, "ACTION")
        hold = et.SubElement(action_hold, "HOLD")
        hold.set("HoldUntilDate", today)

    with sub_xml.open ("wb") as submission_file:
        dom = minidom.parseString(et.tostring(submission, encoding="utf-8"))
        submission_file.write(dom.toprettyxml().encode("utf-8"))

    return sub_xml


def generate_genome_manifest(genome_info, study, manifests_root, alias_to_sample, genome_type, tpa):
    manifest_path = manifests_root / f'{genome_info["genome_name"]}.manifest'

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
                f"This is a {tpa_addition}bin derived from the primary whole genome shotgun (WGS) data set "
                f"{genome_info['study']}. This sample represents a {assembly_type} from the metagenomic run{multiple_runs} "
                f"{genome_info['accessions']}."
            ),
        ),
        ("RUN_REF", genome_info["accessions"]),
        ("FASTA", str(Path(genome_info["genome_path"]).resolve()))
    )
    logger.info(f"Writing manifest file (.manifest) for {genome_info['alias']}.")
    with manifest_path.open("w") as outfile:
        for (k, v) in values:
            manifest = f"{k}\t{v}\n"
            outfile.write(manifest)
        if tpa:
            outfile.write("TPA\ttrue\n")


class GenomeUpload:
    def __init__( self, args):
        self.genome_type = "bins" if args["bins"] else "MAGs"
        self.live = args["live"]
        self.private = args["private"]

        self.tpa = args["tpa"]
        self.centre_name = args["centre_name"]
        self.force = args["force"]

        self.work_dir = Path(args["out"])
        self.upload_dir = self.generate_genomes_upload_dir()
        self.upload_study = args["upload_study"]
        self.genome_metadata = Path(args["genome_info"])
        self.test_suffix = args["test_suffix"]

    def validate_metadata_tsv(self):
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
            raise ValueError(f"The following mandatory fields have missing values in the input file: {', '.join(missing_values)}")

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
                f"Run accessions are not correctly formatted for the following genomes: {', '.join(mismatching_accessions.values)}"
            )

        # check whether completeness and contamination are floats
        try:
            pd.to_numeric(metadata["completeness"])
            pd.to_numeric(metadata["contamination"])
            pd.to_numeric(metadata["genome_coverage"])
        except ValueError:
            raise ValueError("Completeness, contamination and coverage values should be formatted as floats")

        # check whether all co-assemblies have more than one run associated and vice versa
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
        invalid_metagenomes = set(metadata['metagenome'][metadata['metagenome'].map(lambda m: m not in METAGENOMES)])
        if invalid_metagenomes:
            logging.info(f"Invalid metagenomes: {', '.join(invalid_metagenomes)}")
            raise ValueError("Metagenomes associated with each genome need to belong to ENA's approved metagenomes list.")

        # do provided file paths exist?
        missing_paths = metadata['genome_path'].map(lambda p: not Path(p).exists())
        missing_path_set = set(metadata.loc[missing_paths, 'genome_path'])
        if missing_path_set:
            logging.info(f"Missing genome paths: {', '.join(missing_path_set)}")
            raise ValueError("Some genome paths do not exist.")

        # create dictionary while checking genome name uniqueness
        uniqueness = metadata["genome_name"].nunique() == metadata["genome_name"].size
        if uniqueness:
            # for test submissions we add a timestamp to allow for more than one submission a day (ENA would block them otherwise)
            if not self.live:
                if self.test_suffix:
                    unique_test_identifier = self.test_suffix
                else:
                    timestamp = str(int(dt.timestamp(dt.now())))
                    unique_test_identifier = timestamp
                unique_names = [row["genome_name"] + "_" + unique_test_identifier for index, row in metadata.iterrows()]
                metadata["unique_genome_name"] = unique_names
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

        backup_file = self.upload_dir / "ENA_backup.json"
        write_mode = "r+"
        if not backup_file.exists() or self.force:
            write_mode = "w"

        counter = 0

        with backup_file.open(write_mode) as file:
            backup_dict = {}
            if not write_mode == "w":
                try:
                    backup_dict = json.load(file)
                    temp_dict = dict(backup_dict)
                    logger.info(f"A backup file {backup_file} for ENA sample metadata has been found.")
                except json.decoder.JSONDecodeError:
                    pass
            for s in study_set:
                ena_query = EnaQuery(s, "study", self.private)
                study_info = ena_query.build_query()
                project_description = study_info["study_description"]
                if not project_description:
                    project_description = study_info["study_title"]

                ena_query = EnaQuery(s, "study_runs", self.private)
                ena_info = ena_query.build_query()
                if not ena_info:
                    raise IOError(f"No runs found on ENA for project {s}.")

                for run, item in enumerate(ena_info):
                    run_accession = ena_info[run]["run_accession"]
                    if run_accession not in backup_dict:
                        if run_accession in runs_set:
                            sample_accession = ena_info[run]["sample_accession"]
                            ena_query = EnaQuery(sample_accession, "sample", self.private)
                            sample_info = ena_query.build_query()

                            location = sample_info["location"]
                            latitude, longitude = "missing: third party data", "missing: third party data"

                            if location:
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

                                if latitude != "missing: third party data":
                                    try:
                                        latitude = "{:.{}f}".format(round(float(latitude), GEOGRAPHY_DIGIT_COORDS), GEOGRAPHY_DIGIT_COORDS)
                                    except ValueError:
                                        raise IOError("Latitude could not be parsed. Check metadata for run {}.".format(run_accession))

                                if longitude != "missing: third party data":
                                    try:
                                        longitude = "{:.{}f}".format(
                                            round(float(longitude), GEOGRAPHY_DIGIT_COORDS), GEOGRAPHY_DIGIT_COORDS
                                        )
                                    except ValueError:
                                        raise IOError("Longitude could not be parsed. Check metadata for run {}.".format(run_accession))

                            country = sample_info["country"].split(":")[0]
                            if country not in GEOGRAPHIC_LOCATIONS:
                                country = "not provided"

                            collection_date = sample_info["collection_date"]
                            if collection_date.lower() in [
                                "not collected",
                                "not provided",
                                "restricted access",
                                "missing: control sample",
                                "missing: sample group",
                                "missing: synthetic construct",
                                "missing: lab stock",
                                "missing: third party data",
                                "missing: data agreement established pre-2023",
                                "missing: endangered species",
                                "missing: human-identifiable",
                            ]:
                                collection_date = collection_date.lower()
                            if not collection_date or collection_date.lower() == "missing" or collection_date.lower() in [
                                "not available", "na"]:
                                collection_date = "missing: third party data"


                            temp_dict[run_accession] = {
                                "instrumentModel": ena_info[run]["instrument_model"],
                                "collectionDate": collection_date,
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
        upload_dir = self.work_dir / upload_name
        upload_dir.mkdir(parents=True, exist_ok=True)
        return upload_dir

    def create_genome_dictionary(self, samples_xml):
        logger.info("Retrieving data for MAG submission...")

        genome_info = self.extract_genomes_info()

        self.extract_ena_info(genome_info)
        logger.info("Writing genome registration XML...")

        write_genomes_xml(genome_info, samples_xml, self.genome_type, self.centre_name, self.tpa)
        logger.info("All files have been written to " + str(self.upload_dir))

        return genome_info

    def genome_upload(self):
        if not self.live:
            logger.warning("Warning: genome submission is not in live mode, " + "files will be validated, but not uploaded.")
        samples_xml = self.upload_dir / "genome_samples.xml"
        submission_xml = self.upload_dir / "submission.xml"
        genome_info, manifest_info = {}, {}

        # submission xml existence
        if not submission_xml.exists()  or self.force:
            submission_xml = write_submission_xml(self.upload_dir, self.centre_name, False)

        # sample xml generation or recovery
        genome_info = self.create_genome_dictionary(samples_xml)

        accessions_gen = "registered_MAGs.tsv"
        if self.genome_type == "bins":
            accessions_gen = accessions_gen.replace("MAG", "bin")
        if not self.live:
            accessions_gen = accessions_gen.replace(".tsv", "_test.tsv")

        accessions_file = self.upload_dir / accessions_gen

        logger.info("Registering genome samples XMLs...")
        ena_submit = EnaSubmit(samples_xml, submission_xml, len(genome_info), self.live)
        alias_accession_map = ena_submit.handle_genomes_registration()

        if len(alias_accession_map) == len(genome_info):
            # all genomes were registered
            save_accessions(alias_accession_map, accessions_file, "w")
        else:
            if len(alias_accession_map) > 0:
                # exclude those from XML
                if samples_xml.exists():
                    samples_xml.unlink(missing_ok=True)
                filtered_genome_info = {k: v for k, v in genome_info.items() if k not in alias_accession_map}
                logger.info("Re-writing genome registration XML...")
                write_genomes_xml(filtered_genome_info, samples_xml, self.genome_type, self.centre_name, self.tpa)
                logger.info("Registering new genome samples XMLs...")
                ena_submit_new = EnaSubmit(samples_xml, submission_xml, len(filtered_genome_info), self.live)
                new_alias_accession_map = ena_submit_new.handle_genomes_registration()
                if len(new_alias_accession_map) == len(filtered_genome_info):
                    # all new genomes were registered
                    save_accessions(new_alias_accession_map, accessions_file, "a")
                    alias_accession_map.update(new_alias_accession_map)
                else:
                    raise Exception("An error occurred during the registration step.")
            else:
                raise Exception("Some genomes could not be submitted to ENA. Please, check the errors above.")

        # manifests creation
        if not self.live:
            manifest_dir = self.upload_dir / "manifests_test"
        else:
            manifest_dir = self.upload_dir / "manifests"
        manifest_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Generating manifest files...")

        manifest_info = compute_manifests(genome_info)

        for m in manifest_info:
            generate_genome_manifest(
                manifest_info[m], self.upload_study, manifest_dir, alias_accession_map, self.genome_type, self.tpa
            )

__version__ = importlib.metadata.version("genome_uploader")

@click.command()
@click.version_option(__version__, message="genome_uploader %(version)s")
@click.option("-u", "--upload_study", required=True, help="Study accession for genomes upload")
@click.option("--genome_info", type=click.Path(exists=True), required=True, help="Genomes metadata file")
@click.option("-m", "--mags", is_flag=True, help="Select for MAG upload")
@click.option("-b", "--bins", is_flag=True, help="Select for bin upload")
@click.option("--out", type=click.Path(), default=Path.cwd(), help="Output folder. Default: working directory")
@click.option("--force", is_flag=True, help="Forces file reset")
@click.option("--live", is_flag=True, help="Uploads on ENA. Omitting this allows validation only")
@click.option(
    "--test-suffix",
    type=str,
    required=False,
    help="Add suffix (for example, date or time) to generate unique submissions to test server. "
    "Use only for test mode (without --live)! Default: script execution timestamp.",
)
@click.option("--tpa", is_flag=True, help="Select if uploading TPA-generated genomes")
@click.option("--centre_name", type=str, help="Centre name for submission")
@click.option(
    "--private",
    is_flag=True,
    required=False,
    default=False,
    help="If data is private",)


def main(upload_study, genome_info, mags, bins, out, force, live, test_suffix, tpa, centre_name, private):
    if mags == bins:
        raise click.UsageError("Must specify only one of --mags or --bins (not both)")

    if test_suffix and live:
        raise click.UsageError("--test-suffix cannot be used with --live")

    args = {
        "upload_study": upload_study,
        "genome_info": genome_info,
        "mags": mags,
        "bins": bins,
        "out": out,
        "force": force,
        "live": live,
        "test_suffix": test_suffix,
        "tpa": tpa,
        "centre_name": centre_name,
        "private": private,
    }

    ena_upload = GenomeUpload(args)

    if not ena_upload.live:
        logger.warning("Warning: genome submission is not in live mode, files will be validated, but not uploaded.")

    ena_upload.genome_upload()


if __name__ == "__main__":
    main()
    logger.info("Completed")
