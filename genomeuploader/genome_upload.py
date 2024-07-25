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

import os
import sys
import logging
import argparse
import re
import json
import pandas as pd
from datetime import date, datetime as dt
from dotenv import load_dotenv
from pathlib import Path

import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import requests

from ena import ENA

from constants import METAGENOMES, GEOGRAPHIC_LOCATIONS, MQ, HQ

logging.basicConfig(level=logging.DEBUG)

logger = logging.getLogger(__name__)

ena = ENA()

GEOGRAPHY_DIGIT_COORDS = 8

'''
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
'''
def read_and_cleanse_metadata_tsv(inputFile, genomeType, live):
    logger.info('Retrieving info for genomes to submit...')

    binMandatoryFields = ["genome_name", "accessions",
        "assembly_software", "binning_software",
        "binning_parameters", "stats_generation_software", "NCBI_lineage",
        "broad_environment", "local_environment", "environmental_medium", "metagenome",
        "co-assembly", "genome_coverage", "genome_path"]
    MAGMandatoryFields = ["rRNA_presence", "completeness", "contamination"]

    allFields = MAGMandatoryFields + binMandatoryFields
    metadata = pd.read_csv(inputFile, sep='\t', usecols=allFields)

    # make sure there are no empty cells
    cleanColumns = list(metadata.dropna(axis=1))
    if genomeType == "MAGs":
        missingValues = [item for item in allFields if item not in cleanColumns]
    else:
        missingValues = [item for item in binMandatoryFields if item not in cleanColumns]

    if missingValues:
        raise ValueError("The following mandatory fields have missing values in " +
            "the input file: {}".format(", ".join(missingValues)))

    # check amount of genomes to register at the same time
    if len(metadata) >= 5000:
        raise ValueError("Genomes need to be registered in batches of 5000 genomes or smaller.")

    # check whether accessions follow the right format
    accessions_regExp = re.compile(r"([E|S|D]R[R|Z]\d{6,})")

    accessionComparison = pd.DataFrame(columns=["genome_name", "attemptive_accessions",
        "correct", "mismatching", "co-assembly"])
    accessionComparison["genome_name"] = metadata["genome_name"]

    accessionComparison["attemptive_accessions"] = metadata["accessions"].map(
        lambda a: len(a.split(',')))

    accessionComparison["correct"] = metadata["accessions"].map(
        lambda a: len(accessions_regExp.findall(a)))

    accessionComparison["mismatching"] = accessionComparison.apply(lambda row:
        True if row["attemptive_accessions"] == row["correct"]
        else None, axis=1).isna()

    mismatchingAccessions = accessionComparison[accessionComparison["mismatching"]]["genome_name"]
    if not mismatchingAccessions.empty:
        raise ValueError("Run accessions are not correctly formatted for the following " +
            "genomes: " + ','.join(mismatchingAccessions.values))

    # check whether completeness and contamination are floats
    try:
        pd.to_numeric(metadata["completeness"])
        pd.to_numeric(metadata["contamination"])
        pd.to_numeric(metadata["genome_coverage"])
    except:
        raise ValueError("Completeness, contamination or coverage values should be formatted as floats")

    # check whether all co-assemblies have more than one run associated and viceversa
    accessionComparison["co-assembly"] = metadata["co-assembly"]
    coassemblyDiscrepancy = metadata[(
        (accessionComparison["correct"] < 2) & (accessionComparison["co-assembly"])) |
        ((accessionComparison["correct"] > 1) & (~accessionComparison["co-assembly"])
        )]["genome_name"]
    if not coassemblyDiscrepancy.empty:
        raise ValueError("The following genomes show discrepancy between number of runs "
            "involved and co-assembly status: " + ','.join(coassemblyDiscrepancy.values))

    # are provided metagenomes part of the accepted metagenome list?
    if False in metadata.apply(lambda row:
        True if row["metagenome"] in METAGENOMES
        else False, axis=1).unique():
        raise ValueError("Metagenomes associated with each genome need to belong to ENA's " +
            "approved metagenomes list.")

    # do provided file paths exist?
    if False in metadata.apply(lambda row:
        True if os.path.exists(row["genome_path"])
        else False, axis =1).unique():
        raise FileNotFoundError("Some genome paths do not exist.")

    # create dictionary while checking genome name uniqueness
    uniqueness = metadata["genome_name"].nunique() == metadata["genome_name"].size
    if uniqueness:
        if not live:
            timestamp = str(int(dt.timestamp(dt.now())))
            timestamp_names = [row["genome_name"] + '_' + timestamp for index, row in metadata.iterrows()]
            metadata["unique_genome_name"] = timestamp_names
            genomeInfo = metadata.set_index("unique_genome_name").transpose().to_dict()
        else:
            genomeInfo = metadata.set_index("genome_name").transpose().to_dict()
    else:
        raise ValueError("Duplicate names found in genome names")

    return genomeInfo

def round_stats(stats):
    newStat = round(float(stats), 2)
    if newStat == 100.0:
        newStat = 100
    if newStat == 0:
        newStat = 0.0

    return newStat

def compute_MAG_quality(completeness, contamination, RNApresence):
    RNApresent = str(RNApresence).lower() in ["true", "yes", "y"]
    quality = MQ
    if float(completeness) >= 90 and float(contamination) <= 5 and RNApresent:
        quality = HQ

    return quality, completeness, contamination

def extract_tax_info(taxInfo):
    # if unclassified, block the execution
    lineage, position, digitAnnotation = taxInfo.split(';'), 0, False
    lineageFirst = lineage[0]
    if "Unclassified " in lineageFirst:
        if "Archaea" in lineageFirst:
            scientificName = "uncultured archaeon"
        elif "Bacteria" in lineageFirst:
            scientificName = "uncultured bacterium"
        elif "Eukaryota" in lineageFirst:
            scientificName = "uncultured eukaryote"
        submittable, taxid, rank = ena.query_scientific_name(scientificName, searchRank=True)
        return taxid, scientificName

    kingdoms = ["Archaea", "Bacteria", "Eukaryota"]
    kingdomTaxa = ["2157", "2", "2759"]

    selectedKingdom, finalKingdom = kingdoms, ""
    if lineage[1].isdigit():
        selectedKingdom = kingdomTaxa
        position = 2
        digitAnnotation = True
    for index, k in enumerate(selectedKingdom):
        if digitAnnotation:
            if k == lineage[position]:
                finalKingdom = selectedKingdom[index]
                break
        else:
            if k in lineage[position]:
                finalKingdom = selectedKingdom[index]
                break

    iterator = len(lineage)-1
    submittable = False
    rank = ""
    while iterator != -1 and not submittable:
        scientificName = lineage[iterator].strip()
        if digitAnnotation:
            if not '*' in scientificName:
                scientificName = ena.query_taxid(scientificName)
            else:
                iterator -= 1
                continue
        elif "__" in scientificName:
            scientificName = scientificName.split("__")[1]
        else:
            raise ValueError("Unrecognised taxonomy format: " + scientificName)
        submittable, taxid, rank = ena.query_scientific_name(scientificName, searchRank=True)

        if not submittable:
            if finalKingdom == "Archaea" or finalKingdom == "2157":
                submittable, scientificName, taxid = extract_Archaea_info(scientificName, rank)
            elif finalKingdom == "Bacteria" or finalKingdom == "2":
                submittable, scientificName, taxid = extract_Bacteria_info(scientificName, rank)
            elif finalKingdom == "Eukaryota" or finalKingdom == "2759":
                submittable, scientificName, taxid = extract_Eukaryota_info(scientificName, rank)
        iterator -= 1

    return taxid, scientificName

def extract_Eukaryota_info(name, rank):
    nonSubmittable = (False, "", 0)

    # Asterisks in given taxonomy suggest the classification might be not confident enough.
    if '*' in name:
        return nonSubmittable

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
                name = name.replace(" sp.", '')
                submittable, taxid = ena.query_scientific_name(name)
                if submittable:
                    return submittable, name, taxid
                else:
                    return nonSubmittable

def extract_Bacteria_info(name, rank):
    if rank == "species":
        name = name
    elif rank == "superkingdom":
        name = "uncultured bacterium".format(name)
    elif rank in ["family", "order", "class", "phylum"]:
        name = "uncultured {} bacterium".format(name)
    elif rank == "genus":
        name = "uncultured {} sp.".format(name)

    submittable, taxid, rank = ena.query_scientific_name(name, searchRank=True)
    if not submittable:
        if rank in ["species", "genus"] and name.lower().endswith("bacteria"):
            name = "uncultured {}".format(name.lower().replace("bacteria", "bacterium"))
        elif rank == "family":
            if name.lower() == "deltaproteobacteria":
                name = "uncultured delta proteobacterium"
        submittable, taxid = ena.query_scientific_name(name)

    return submittable, name, taxid

def extract_Archaea_info(name, rank):
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

    submittable, taxid, rank = ena.query_scientific_name(name, searchRank=True)
    if not submittable:
        if "Candidatus" in name:
            if rank == "phylum":
                name = name.replace("Candidatus ", '')
            elif rank == "family":
                name = name.replace("uncultured ", '')
            submittable, taxid = ena.query_scientific_name(name)

    return submittable, name, taxid

def extract_genomes_info(inputFile, genomeType, live):
    genomeInfo = read_and_cleanse_metadata_tsv(inputFile, genomeType, live)
    for gen in genomeInfo:
        genomeInfo[gen]["accessions"] = genomeInfo[gen]["accessions"].split(',')
        accessionType = "run"
        assembly_regExp = re.compile(r"([E|S|D]RZ\d{6,})")
        if assembly_regExp.findall(genomeInfo[gen]["accessions"][0]):
            accessionType = "assembly"
        genomeInfo[gen]["accessionType"] = accessionType

        genomeInfo[gen]["isolationSource"] = genomeInfo[gen]["metagenome"]

        try:
            (genomeInfo[gen]["MAG_quality"],
            genomeInfo[gen]["completeness"],
            genomeInfo[gen]["contamination"]) = compute_MAG_quality(
                                    str(round_stats(genomeInfo[gen]["completeness"])),
                                    str(round_stats(genomeInfo[gen]["contamination"])),
                                    genomeInfo[gen]["rRNA_presence"])
        except IndexError:
            pass

        if str(genomeInfo[gen]["co-assembly"]).lower() in ["yes", "y", "true"]:
            genomeInfo[gen]["co-assembly"] = True
        else:
            genomeInfo[gen]["co-assembly"] = False

        genomeInfo[gen]["alias"] = gen

        taxID, scientificName = extract_tax_info(genomeInfo[gen]["NCBI_lineage"])
        genomeInfo[gen]["taxID"] = taxID
        genomeInfo[gen]["scientific_name"] = scientificName

    return genomeInfo

def extract_ENA_info(genomeInfo, uploadDir, webin, password):
    logger.info('Retrieving project and run info from ENA (this might take a while)...')

    # retrieving metadata from runs (and runs from assembly accessions if provided)
    allRuns = []
    for g in genomeInfo:
        if genomeInfo[g]["accessionType"] == "assembly":
            derivedRuns = []
            for acc in genomeInfo[g]["accessions"]:
                derivedRuns.append(ena.get_run_from_assembly(acc))
            genomeInfo[g]["accessions"] = derivedRuns
        allRuns.extend(genomeInfo[g]["accessions"])

    runsSet, studySet, samplesDict, tempDict = set(allRuns), set(), {}, {}
    for r in runsSet:
        run_info = ena.get_run(r, webin, password)
        studySet.add(run_info["secondary_study_accession"])
        samplesDict[r] = run_info["sample_accession"]

    if not studySet:
        raise ValueError("No study corresponding to runs found.")

    backupFile = os.path.join(uploadDir, "ENA_backup.json")
    counter = 0
    if not os.path.exists(backupFile):
        with open(backupFile, 'w') as file:
            pass
    with open(backupFile, "r+") as file:
        try:
            backupDict = json.load(file)
            tempDict = dict(backupDict)
            logger.info(f"A backup file {backupFile} for ENA sample metadata has been found.")
        except json.decoder.JSONDecodeError:
            backupDict = {}
        for s in studySet:
            studyInfo = ena.get_study(webin, password, s)
            projectDescription = studyInfo["study_description"]
            if not projectDescription:
                projectDescription = studyInfo["study_title"]

            ENA_info = ena.get_study_runs(s, webin, password)
            if ENA_info == []:
                raise IOError("No runs found on ENA for project {}.".format(s))

            for run, item in enumerate(ENA_info):
                runAccession = ENA_info[run]["run_accession"]
                if runAccession not in backupDict:
                    if runAccession in runsSet:
                        sampleAccession = ENA_info[run]["sample_accession"]
                        sampleInfo = ena.get_sample(sampleAccession, webin, password)

                        location = sampleInfo["location"]
                        latitude, longitude = "not provided", "not provided"
                        if 'N' in location:
                            latitude = location.split('N')[0].strip()
                            longitude = location.split('N')[1].strip()
                        elif 'S' in location:
                            latitude = '-' + location.split('S')[0].strip()
                            longitude = location.split('S')[1].strip()

                        if 'W' in longitude:
                            longitude = '-' + longitude.split('W')[0].strip()
                        elif longitude.endswith('E'):
                            longitude = longitude.split('E')[0].strip()

                        if latitude != "not provided":
                            try:
                                latitude = "{:.{}f}".format(round(float(latitude), GEOGRAPHY_DIGIT_COORDS), GEOGRAPHY_DIGIT_COORDS)
                            except ValueError:
                                raise IOError("Latitude could not be parsed. Check metadata for run {}.".format(runAccession))

                        if longitude != "not provided":
                            try:
                                longitude = "{:.{}f}".format(round(float(longitude), GEOGRAPHY_DIGIT_COORDS), GEOGRAPHY_DIGIT_COORDS)
                            except ValueError:
                                raise IOError("Longitude could not be parsed. Check metadata for run {}.".format(runAccession))

                        country = sampleInfo["country"].split(':')[0]
                        if not country in GEOGRAPHIC_LOCATIONS:
                            country = "not provided"

                        collectionDate = sampleInfo["collection_date"]
                        if not collectionDate or collectionDate.lower() == "missing" or collectionDate.lower() in \
                                ["not applicable", "not available", "na"]:
                            collectionDate = "not provided"

                        tempDict[runAccession] = {
                            "instrumentModel" : ENA_info[run]["instrument_model"],
                            "collectionDate" : collectionDate,
                            "country" : country,
                            "latitude" : latitude,
                            "longitude" : longitude,
                            "projectDescription" : projectDescription,
                            "study" : s,
                            "sampleAccession" : samplesDict[runAccession]
                        }
                        counter += 1

                        if (counter%10 == 0) or (len(runsSet) - len(backupDict) == counter):
                            file.seek(0)
                            file.write(json.dumps(tempDict))
                            file.truncate()
    tempDict = {**tempDict, **backupDict}
    combine_ENA_info(genomeInfo, tempDict)

def multipleElementSet(metadataList):
    return len(set(metadataList))>1

def combine_ENA_info(genomeInfo, ENADict):
    for g in genomeInfo:
        # TODO: optimise all the part below
        if genomeInfo[g]["co-assembly"]:
            instrumentList, collectionList, countryList = [], [], []
            studyList, descriptionList, samplesList = [], [], []
            longList, latitList = [], []
            for run in genomeInfo[g]["accessions"]:
                instrumentList.append(ENADict[run]["instrumentModel"])
                collectionList.append(ENADict[run]["collectionDate"])
                countryList.append(ENADict[run]["country"])
                studyList.append(ENADict[run]["study"])
                descriptionList.append(ENADict[run]["projectDescription"])
                samplesList.append(ENADict[run]["sampleAccession"])
                longList.append(ENADict[run]["longitude"])
                latitList.append(ENADict[run]["latitude"])

            genomeInfo[g]["study"] = studyList[0]
            genomeInfo[g]["description"] = descriptionList[0]

            instrument = instrumentList[0]
            if multipleElementSet(instrumentList):
                instrument = ','.join(instrumentList)
            genomeInfo[g]["sequencingMethod"] = instrument

            collectionDate = collectionList[0]
            if multipleElementSet(collectionList):
                collectionDate = "not provided"
            if collectionDate.lower() in ["not applicable", "not available", "na"]:
                collectionDate = "not provided"
            genomeInfo[g]["collectionDate"] = collectionDate

            country = countryList[0]
            if multipleElementSet(countryList):
                country = "not applicable"
            genomeInfo[g]["country"] = country

            latitude = latitList[0]
            if multipleElementSet(latitList):
                latitude = "not provided"
            try:
                genomeInfo[g]["latitude"] = str(round(float(latitude), GEOGRAPHY_DIGIT_COORDS))
            except ValueError:
                genomeInfo[g]["latitude"] = "not provided"

            longitude = longList[0]
            if multipleElementSet(longList):
                longitude = "not provided"
            try:
                genomeInfo[g]["longitude"] = str(round(float(longitude), GEOGRAPHY_DIGIT_COORDS))
            except ValueError:
                genomeInfo[g]["longitude"] = "not provided"

            samples = samplesList[0]
            if multipleElementSet(samplesList):
                samples = ','.join(samplesList)
            genomeInfo[g]["sample_accessions"] = samples
        else:
            run = genomeInfo[g]["accessions"][0]
            genomeInfo[g]["sequencingMethod"] = ENADict[run]["instrumentModel"]
            if ENADict[run]["collectionDate"].lower() in ["not applicable", "not available", "na"]:
                genomeInfo[g]["collectionDate"] = "not provided"
            else:
                genomeInfo[g]["collectionDate"] = ENADict[run]["collectionDate"]
            genomeInfo[g]["study"] = ENADict[run]["study"]
            genomeInfo[g]["description"] = ENADict[run]["projectDescription"]
            genomeInfo[g]["sample_accessions"] = ENADict[run]["sampleAccession"]
            genomeInfo[g]["country"] = ENADict[run]["country"]
            genomeInfo[g]["longitude"] = ENADict[run]["longitude"]
            genomeInfo[g]["latitude"] = ENADict[run]["latitude"]

        genomeInfo[g]["accessions"] = ','.join(genomeInfo[g]["accessions"])


def getAccessions(accessionsFile):
    accessionDict = {}
    with open(accessionsFile, 'r') as f:
        for line in f:
            line = line.split('\t')
            alias = line[0]
            accession = line[1].rstrip('\n')
            accessionDict[alias] = accession

    return accessionDict

def saveAccessions(aliasAccessionDict, accessionsFile, writeMode):
    with open(accessionsFile, writeMode) as f:
        for elem in aliasAccessionDict:
            f.write("{}\t{}\n".format(elem, aliasAccessionDict[elem]))

def create_manifest_dictionary(run, alias, assemblySoftware, sequencingMethod,
    MAGpath, gen, study, coverage, isCoassembly):
    manifestDict = {
        "accessions" : run,
        "alias" : alias,
        "assembler" : assemblySoftware,
        "sequencingMethod" : sequencingMethod,
        "genome_path" : MAGpath,
        "genome_name" : gen,
        "study" : study,
        "coverageDepth" : coverage,
        "co-assembly" : isCoassembly
    }

    return manifestDict

def compute_manifests(genomes):
    manifestInfo = {}
    for g in genomes:
        manifestInfo[g] = create_manifest_dictionary(genomes[g]["accessions"],
            genomes[g]["alias"], genomes[g]["assembly_software"],
            genomes[g]["sequencingMethod"], genomes[g]["genome_path"], g,
            genomes[g]["study"], genomes[g]["genome_coverage"], genomes[g]["co-assembly"])

    return manifestInfo

def get_study_from_xml(sample):
    description = sample.childNodes[5].childNodes[0].data
    study = description.split(' ')[-1][:-1]

    return study

def recover_info_from_xml(genomeDict, sample_xml, live_mode):
    logger.info("Retrieving data for genome submission...")

    # extract list of genomes (samples) to be registered
    xml_structure = minidom.parse(sample_xml)
    samples = xml_structure.getElementsByTagName("SAMPLE")

    for s in samples:
        study = get_study_from_xml(s)

        # extract alias from xml and find a match with genomes the user is uploading
        XMLalias = s.attributes["alias"].value
        if not live_mode:         # remove time stamp if test mode is selected
            aliasSplit = XMLalias.split("_")
            XMLalias = '_'.join(aliasSplit[:-1])
        for gen in genomeDict:
            # if match is found, associate attributes listed in the xml file
            # with genomes to upload
            if XMLalias == gen:
                if not live_mode:
                    currentTimestamp = str(int(dt.timestamp(dt.now())))
                    XMLalias = gen + '_' + currentTimestamp
                    s.attributes["alias"].value = XMLalias
                    sampleTitle = s.getElementsByTagName("TITLE")[0]
                    sampleTitleValue = sampleTitle.firstChild.nodeValue.split("_")
                    sampleTitleValue[-1] = currentTimestamp
                    newSampleTitle = '_'.join(sampleTitleValue)
                    s.getElementsByTagName("TITLE")[0].firstChild.replaceWholeText(newSampleTitle)
                attributes = s.childNodes[7].getElementsByTagName("SAMPLE_ATTRIBUTE")
                seqMethod, assSoftware = "", ""
                for a in attributes:
                    tagElem = a.getElementsByTagName("TAG")
                    tag = tagElem[0].childNodes[0].nodeValue
                    if tag == "sequencing method":
                        seqMethodElem = a.getElementsByTagName("VALUE")
                        seqMethod = seqMethodElem[0].childNodes[0].nodeValue
                    elif tag == "assembly software":
                        assSoftwareElem = a.getElementsByTagName("VALUE")
                        assSoftware = assSoftwareElem[0].childNodes[0].nodeValue
                    if not seqMethod == "" and not assSoftware == "":
                        break

                genomeDict[gen]["accessions"] = ','.join(genomeDict[gen]["accessions"])
                genomeDict[gen]["alias"] = XMLalias
                genomeDict[gen]["assembly_software"] = assSoftware
                genomeDict[gen]["sequencingMethod"] = seqMethod
                genomeDict[gen]["study"] = study
                break

    if not live_mode:
        for s in samples:
            xml_structure.firstChild.appendChild(s)

        with open(sample_xml, 'wb') as f:
            dom_string = xml_structure.toprettyxml().encode("utf-8")
            dom_string = b'\n'.join([s for s in dom_string.splitlines() if s.strip()])
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
    ET.SubElement(new_sample_attr, 'TAG').text = tag
    ET.SubElement(new_sample_attr, 'VALUE').text = value

    if units:
        ET.SubElement(new_sample_attr, 'UNITS').text = units

def write_genomes_xml(genomes, xml_path, genomeType, centreName, tpa):
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

    checklist, assemblyType = "ERC000047", "Metagenome-assembled genome"
    if genomeType == "bins":
        checklist = "ERC000050"
        assemblyType = "binned metagenome"

    constant_sample_attributes = [
        # tag - value
        ["taxonomic identity marker", "multi-marker approach"],
        ["investigation type", "metagenome-assembled genome"],
        ["ENA-CHECKLIST", checklist],
    ]

    tpaDescription = ""
    if tpa:
        tpaDescription = "Third Party Annotation (TPA) "

    sample_set = ET.Element("SAMPLE_SET")

    for g in genomes:
        plural = ""
        if genomes[g]["co-assembly"]:
            plural = 's'
        description = ("This sample represents a {}{} assembled from the "
            "metagenomic run{} {} of study {}.".format(tpaDescription,
            assemblyType, plural, genomes[g]["accessions"],
            genomes[g]["study"]))

        sample = ET.SubElement(sample_set, "SAMPLE")
        sample.set("alias", genomes[g]["alias"])
        sample.set("center_name", centreName)

        ET.SubElement(sample, 'TITLE').text = ("{}: {}".format(assemblyType,
            genomes[g]["alias"]))
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

    with open(xml_path, 'wb') as f:
        dom = minidom.parseString(
            ET.tostring(sample_set, encoding="utf-8")
        )
        f.write(dom.toprettyxml().encode("utf-8"))

def write_submission_xml(upload_dir, centre_name, study=True):
    today = str(date.today())
    sub_xml = os.path.join(upload_dir, 'submission.xml')

    submission = ET.Element('SUBMISSION')
    submission.set('center_name', centre_name)

    # template
    actions = ET.SubElement(submission, 'ACTIONS')
    action_sub = ET.SubElement(actions, 'ACTION')
    ET.SubElement(action_sub, 'ADD')

    # attributes: function and hold date
    if study:
        action_hold = ET.SubElement(actions, 'ACTION')
        hold = ET.SubElement(action_hold, 'HOLD')
        hold.set('HoldUntilDate', today)

    with open(sub_xml, 'wb') as submission_file:
        dom = minidom.parseString(
            ET.tostring(submission, encoding="utf-8")
        )
        submission_file.write(dom.toprettyxml().encode("utf-8"))

    return sub_xml

def generate_genome_manifest(genomeInfo, study, manifestsRoot, aliasToSample, genomeType, tpa):
    manifest_path = os.path.join(manifestsRoot, f'{genomeInfo["genome_name"]}.manifest')

    tpaAddition, multipleRuns = "", ""
    if tpa:
        tpaAddition = "Third Party Annotation (TPA) "
    if genomeInfo["co-assembly"]:
        multipleRuns = "s"
    assemblyType = "Metagenome-Assembled Genome (MAG)"
    if genomeType == "bins":
        assemblyType = "binned metagenome"

    values = (
        ('STUDY', study),
        ('SAMPLE', aliasToSample[genomeInfo["alias"]]),
        ('ASSEMBLYNAME', genomeInfo["alias"]),
        ('ASSEMBLY_TYPE', assemblyType),
        ('COVERAGE', genomeInfo["coverageDepth"]),
        ('PROGRAM', genomeInfo["assembler"]),
        ('PLATFORM', genomeInfo["sequencingMethod"]),
        ('MOLECULETYPE', "genomic DNA"),
        ('DESCRIPTION', ("This is a {}bin derived from the primary whole genome "
            "shotgun (WGS) data set {}. This sample represents a {} from the "
            "metagenomic run{} {}.".format(tpaAddition, genomeInfo["study"],
            assemblyType, multipleRuns, genomeInfo["accessions"]))),
        ('RUN_REF', genomeInfo["accessions"]),
        ('FASTA', os.path.abspath(genomeInfo["genome_path"]))
    )
    logger.info("Writing manifest file (.manifest) for {}.".format(genomeInfo["alias"]))
    with open(manifest_path, "w") as outfile:
        for (k, v) in values:
            manifest = f'{k}\t{v}\n'
            outfile.write(manifest)
        if tpa:
            outfile.write("TPA\ttrue\n")

def main():
    ENA_uploader = GenomeUpload()

    if not ENA_uploader.live:
        logger.warning("Warning: genome submission is not in live mode, " +
            "files will be validated, but not uploaded.")

    xmlGenomeFile, xmlSubFile = "genome_samples.xml", "submission.xml"
    samples_xml = os.path.join(ENA_uploader.upload_dir, xmlGenomeFile)
    submissionXmlPath = os.path.join(ENA_uploader.upload_dir, xmlSubFile)
    submission_xml = submissionXmlPath
    genomes, manifestInfo = {}, {}

    # submission xml existence
    if not os.path.exists(submissionXmlPath):
        submission_xml = write_submission_xml(ENA_uploader.upload_dir, ENA_uploader.centre_name, False)

    # sample xml generation or recovery
    genomes = ENA_uploader.create_genome_dictionary(samples_xml)

    # manifests creation
    manifestDir = os.path.join(ENA_uploader.upload_dir, "manifests")
    os.makedirs(manifestDir, exist_ok=True)

    accessionsgen = "registered_MAGs.tsv"
    if ENA_uploader.genomeType == "bins":
        accessionsgen = accessionsgen.replace("MAG", "bin")
    if not ENA_uploader.live:
        accessionsgen = accessionsgen.replace(".tsv", "_test.tsv")

    accessionsFile = os.path.join(ENA_uploader.upload_dir, accessionsgen)
    save = False
    writeMode = 'a'
    if os.path.exists(accessionsFile):
        if not ENA_uploader.live:
            save = True
            if ENA_uploader.force:
                writeMode = 'w'
        if not save:
            logger.info("Genome samples already registered, reading ERS accessions...")
            aliasToNewSampleAccession = getAccessions(accessionsFile)
    else:
        save = True

    if save:
        logger.info("Registering genome samples XMLs...")
        aliasToNewSampleAccession = ena.handle_genomes_registration(samples_xml,
            submission_xml, ENA_uploader.username, ENA_uploader.password, ENA_uploader.live)
        saveAccessions(aliasToNewSampleAccession, accessionsFile, writeMode)

    logger.info("Generating manifest files...")

    manifestInfo = compute_manifests(genomes)

    for m in manifestInfo:
        generate_genome_manifest(manifestInfo[m], ENA_uploader.upStudy,
            manifestDir, aliasToNewSampleAccession, ENA_uploader.genomeType, ENA_uploader.tpa)

class GenomeUpload:
    def __init__(self, argv=sys.argv[1:]):
        self.args = self.parse_args(argv)
        self.upStudy = self.args.upload_study
        self.genomeMetadata = self.args.genome_info
        self.genomeType = "bins" if self.args.bins else "MAGs"
        self.live = True if self.args.live else False

        if self.args.webin and self.args.password:
            self.username = self.args.webin
            self.password = self.args.password
        else:
            # Config file
            user_config = Path.home() / ".genome_uploader.config.env"
            if user_config.exists():
                logger.debug("Loading the env variables from ".format(user_config))
                load_dotenv(str(user_config))
            else:
                cwd_config = Path.cwd() / ".genome_uploader.config.env"
                if cwd_config.exists():
                    logger.debug("Loading the variables from the current directory.")
                    load_dotenv(str(cwd_config))
                else:
                    logger.debug("Trying to load env variables from the .env file")
                    # from a local .env file
                    load_dotenv()

            self.username = os.getenv("ENA_WEBIN")
            self.password = os.getenv("ENA_WEBIN_PASSWORD")

        if not self.username or not self.password:
            logger.error("ENA Webin username or password are empty")
            sys.exit(1)

        self.tpa = True if self.args.tpa else False
        self.centre_name = self.args.centre_name
        self.force = True if self.args.force else False

        workDir = self.args.out if self.args.out else os.getcwd()
        self.upload_dir = self.generate_genomes_upload_dir(workDir, self.genomeType)

    def parse_args(self, argv):
        parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description="Create xmls and manifest files for genome upload to ENA")

        parser.add_argument('-u', '--upload_study', type=str, help="Study accession for genomes upload")
        parser.add_argument('--genome_info', type=str, required=True, help="Genomes metadata file")

        genomeType = parser.add_mutually_exclusive_group(required=True)
        genomeType.add_argument('-m', '--mags', action='store_true', help="Select for MAG upload")
        genomeType.add_argument('-b', '--bins', action='store_true', help="Select for bin upload")

        parser.add_argument('--out', type=str, help="Output folder. Default: working directory")
        parser.add_argument('--force', action='store_true', help="Forces reset of sample xml's backups")
        parser.add_argument('--live', action='store_true', help="Uploads on ENA. Omitting this " +
            "option allows to validate samples beforehand")
        parser.add_argument('--tpa', action='store_true', help="Select if uploading TPA-generated genomes")

        # Users can provide their credentials and centre name manually or using a config file
        parser.add_argument('--webin', required=False, help="Webin id")
        parser.add_argument('--password', required=False, help="Webin password")
        parser.add_argument('--centre_name', required=False, help="Name of the centre uploading genomes")

        args = parser.parse_args(argv)

        if not args.upload_study:
            raise ValueError("No project selected for genome upload [-u, --upload_study].")

        if not os.path.exists(args.genome_info):
            raise FileNotFoundError('Genome metadata file "{}" does not exist'.format(args.genome_info))

        return args

    def generate_genomes_upload_dir(self, dir, genomeType):
        uploadName = "MAG_upload"
        if genomeType == "bins":
            uploadName = uploadName.replace("MAG", "bin")
        upload_dir = os.path.join(dir, uploadName)
        os.makedirs(upload_dir, exist_ok=True)
        return upload_dir

    def create_genome_dictionary(self, samples_xml):
        logger.info('Retrieving data for MAG submission...')

        genomeInfo = extract_genomes_info(self.genomeMetadata, self.genomeType, self.live)

        if not os.path.exists(samples_xml) or self.force:
            extract_ENA_info(genomeInfo, self.upload_dir, self.username, self.password)
            logger.info("Writing genome registration XML...")

            write_genomes_xml(genomeInfo, samples_xml, self.genomeType,
                              self.centre_name, self.tpa)
            logger.info("All files have been written to " + self.upload_dir)
        else:
            recover_info_from_xml(genomeInfo, samples_xml, self.live)

        return genomeInfo

if __name__ == "__main__":
    main()
    logger.info('Completed')
