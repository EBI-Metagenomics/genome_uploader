#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2017-2025 EMBL - European Bioinformatics Institute
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
import requests
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TaxonFinder:
    def __init__(self, start_lineage):
        self.lineage = start_lineage
        self.taxid, self.scientific_name = self.extract_tax_info(self.lineage)

    def query_taxid(self, taxid):
        """
        Queries ENA taxonomy API for a scientific name given a taxid.
        Args:
            taxid (str or int): NCBI taxonomic ID.
        Returns:
            str: Scientific name corresponding to the taxid, or empty string if not found.
        """
        url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{taxid}"
        response = requests.get(url)

        try:
            # Will raise exception if response status code is non-200
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logger.error(f"Request failed {url} with error {e}")
            return False
        res = response.json()
        return res.get("scientificName", "")

    def query_scientific_name(self, scientific_name, search_rank=False):
        """
        Queries ENA taxonomy API for a scientific name and optionally its rank.
        If the name exists, it checks if it's submittable and retrieves its taxid.
        Args:
            scientific_name (str): Scientific name to query.
            search_rank (bool): If True, also return taxonomic rank.
        Returns:
            tuple or bool: If search_rank is True, returns (submittable, taxid,
            rank). Otherwise, returns (submittable, taxid).
        """
        url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{scientific_name}"
        response = requests.get(url)

        try:
            # Will raise exception if response status code is non-200
            response.raise_for_status()
        except requests.exceptions.HTTPError:
            if search_rank:
                return False, "", ""
            else:
                return False, ""

        try:
            res = response.json()[0]
        except IndexError:
            if search_rank:
                return False, "", ""
            else:
                return False, ""

        submittable = res.get("submittable", "").lower() == "true"
        taxid = res.get("taxId", "")
        rank = res.get("rank", "")

        if search_rank:
            return submittable, taxid, rank
        else:
            return submittable, taxid

    def extract_tax_info(self, tax_info: str) -> tuple[int, str]:
        """
        Extracts taxonomic information from a lineage string and determines
        a valid scientific name and taxid for ENA submission. This function
        parses a semicolon-separated lineage string (with taxon names or ids),
        identifies the kingdom, and iteratively searches for a submittable
        scientific name and taxid using ENA queries. It handles special
        cases for unclassified lineages and rolls up the taxonomy tree if
        needed, applying custom rules for Archaea, Bacteria, and Eukaryota.
        Args:
            tax_info (str): NCBI taxonomic lineage string (with semicolons
            separating taxon names or ids i.e. list of strings or integers).
        Returns:
            tuple[int, str]: (taxid, scientific_name)
        """
        taxid = ''
        lineage = tax_info.split(";")

        # if unclassified, block the execution and get the official name for the kingdom
        lineage_first = lineage[0]
        scientific_name = ''
        if "Unclassified " in lineage_first:
            if "Archaea" in lineage_first:
                scientific_name = "uncultured archaeon"
            elif "Bacteria" in lineage_first:
                scientific_name = "uncultured bacterium"
            elif "Eukaryota" in lineage_first:
                scientific_name = "uncultured eukaryote"
            submittable, taxid, rank = self.query_scientific_name(scientific_name, search_rank=True)
            return taxid, scientific_name

        kingdoms = ["Archaea", "Bacteria", "Eukaryota"]
        kingdom_taxa = ["2157", "2", "2759"]

        # If the provided taxonomic annotation is a string, the first element
        # in the list represents the kingdom. Otherwise, if it was provided as a
        # list of integers, the first element of the lineage will be the NCBI
        # root id, and the second element will indicate the kingdom. Once the
        # kingdom slot is identified in the lineage, the script then checks its
        # value against the existing kingdom list.
        digit_annotation = False
        kingdom_position_lineage = 0
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

        # start iterating from the most specific taxonomic level and roll up a
        # level in the taxonomy tree if the identified name is unsubmittable
        iterator = len(lineage) - 1
        submittable = False
        while iterator != -1 and not submittable:
            scientific_name = lineage[iterator].strip()
            if digit_annotation:
                if "*" not in scientific_name:
                    scientific_name = self.query_taxid(scientific_name)
                else:
                    iterator -= 1
                    continue
            # format if using format like "k__;...g__;s__"
            elif "__" in scientific_name:
                scientific_name = scientific_name.split("__")[1]
            else:
                raise ValueError("Unrecognised taxonomy format: " + scientific_name)
            submittable, taxid, rank = self.query_scientific_name(scientific_name, search_rank=True)

            if not submittable:
                if final_kingdom == "Archaea" or final_kingdom == "2157":
                    submittable, scientific_name, taxid = self.extract_archaea_info(scientific_name, rank)
                elif final_kingdom == "Bacteria" or final_kingdom == "2":
                    submittable, scientific_name, taxid = self.extract_bacteria_info(scientific_name, rank)
                elif final_kingdom == "Eukaryota" or final_kingdom == "2759":
                    submittable, scientific_name, taxid = self.extract_eukaryota_info(scientific_name, rank)
            iterator -= 1

        return taxid, scientific_name

    def extract_eukaryota_info(self, name: str, rank: str) -> tuple[bool, str, int]:
        """
        Extracts and formats eukaryotic taxonomic information for ENA submission.
        It tries querying ENA API to verify whether custom names and/or exceptions
        already exist in the NCBI taxonomy for the queried organism without having to
        roll up a level in the taxonomy tree. If the name exists, it labels the
        organism as submittable.
        Args:
            name (str): Taxonomic name.
            rank (str): Taxonomic rank.
        Returns:
            tuple[bool, str, int]: (submittable, scientific_name, taxid)
        """
        non_submittable = (False, "", 0)

        # Asterisks in given taxonomy suggest the classification might be not confident enough.
        if "*" in name:
            return non_submittable

        if rank == "super kingdom":
            name = "uncultured eukaryote"
            submittable, taxid = self.query_scientific_name(name)
            return submittable, name, taxid
        else:
            name = name.capitalize() + " sp."
            submittable, taxid = self.query_scientific_name(name)
            if submittable:
                return submittable, name, taxid
            else:
                name = "uncultured " + name
                submittable, taxid = self.query_scientific_name(name)
                if submittable:
                    return submittable, name, taxid
                else:
                    name = name.replace(" sp.", "")
                    submittable, taxid = self.query_scientific_name(name)
                    if submittable:
                        return submittable, name, taxid
                    else:
                        return non_submittable

    def extract_bacteria_info(self, name: str, rank: str) -> tuple[bool, str, int]:
        """
        Extracts and formats bacterial taxonomic information for ENA submission.
        It tries querying ENA API to verify whether custom names and/or exceptions
        already exist in the NCBI taxonomy for the queried organism without having to
        roll up a level in the taxonomy tree. If the name exists, it labels the
        organism as submittable.
        Args:
            name (str): Taxonomic name.
            rank (str): Taxonomic rank.
        Returns:
            tuple[bool, str, int]: (submittable, scientific_name, taxid)
        """
        if rank == "species":
            name = name
        elif rank == "domain":
            name = "uncultured bacterium"
        elif rank in ["family", "order", "class", "phylum"]:
            name = f"{name} bacterium"
        elif rank == "genus":
            name = f"{name} sp."

        submittable, taxid, rank = self.query_scientific_name(name, search_rank=True)
        if not submittable:
            name = name.lower()
            if rank in ["species", "genus", "family"]:
                if name.endswith("bacteria"):
                    name = f"uncultured {name[:-9]} bacterium" # replace the last instance of 'bacteria' with 'bacterium'
                else:
                    name = f"uncultured {name}"
            elif "alphaproteobacteria" in name:
                name = "uncultured alpha proteobacterium"
            elif "betaproteobacteria" in name:
                name = "uncultured beta proteobacterium"
            elif "gammaproteobacteria" in name:
                name = "uncultured gammaproteobacteria bacterium"
            elif "deltaproteobacteria" in name:
                name = "uncultured delta proteobacterium"
            submittable, taxid = self.query_scientific_name(name)

        return submittable, name, taxid

    def extract_archaea_info(self, name: str, rank: str) -> tuple[bool, str, int]:
        """
        Extracts and formats archaeal taxonomic information for ENA submission.
        It tries querying ENA API to verify whether custom names and/or exceptions
        already exist in the NCBI taxonomy for the queried organism without having to
        roll up a level in the taxonomy tree. If the name exists, it labels the
        organism as submittable.
        Args:
            name (str): Taxonomic name.
            rank (str): Taxonomic rank.
        Returns:
            tuple[bool, str, int]: (submittable, scientific_name, taxid)
        """
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

        submittable, taxid, rank = self.query_scientific_name(name, search_rank=True)
        if not submittable:
            if "Candidatus" in name:
                if rank == "phylum":
                    name = name.replace("Candidatus ", "")
                elif rank == "family":
                    name = name.replace("uncultured ", "")
                submittable, taxid = self.query_scientific_name(name)

        return submittable, name, taxid


def main():
    parser = argparse.ArgumentParser(
        description="Find taxID and scientific name from a lineage string"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file with lineages one per line"
    )
    args = parser.parse_args()
    with open(args.input, 'r') as file_in:
        for line in file_in:
            line = line.strip()
            finder = TaxonFinder(line)
            print(f"Lineage: {line}")
            print(f"Scientific name: {finder.scientific_name}")
            print(f"TaxID: {finder.taxid}")


if __name__ == "__main__":
    main()