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

"""
Everything specific to submitting single-contig genomes as chromosomes: parsing the
optional --single-contig-info TSV, validating its values, writing the ENA chromosome
list file, and a dedicated log for problems found along the way.
"""

import gzip
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from genomeuploader.constants import (
    GENOME_NAME_FIELD,
    CHROMOSOME_NAME_FIELD,
    CHROMOSOME_NAME_REGEX,
    CHROMOSOME_TYPE_FIELD,
    CHROMOSOME_TYPES,
    CHROMOSOME_TOPOLOGY_FIELD,
    CHROMOSOME_TOPOLOGIES,
    CHROMOSOME_LOCATION_FIELD,
    CHROMOSOME_LOCATIONS_LIST,
)


def normalise_na(value) -> str:
    """
    Normalises a metadata cell to a plain string.
    pandas represents empty cells as NaN (a float), which is truthy and breaks
    membership/format checks; this converts NaN/None to an empty string and
    strips surrounding whitespace from everything else.
    Args:
        value: Raw value read from the metadata table.
    Returns:
        str: Empty string for NaN/None, otherwise the stripped string form.
    """
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def load_single_contig_metadata(single_contig_info: Optional[Path]) -> dict:
    """
    Loads and validates the optional single-contig metadata TSV (--single-contig-info).

    Every genome_name listed in this file is submitted as a single-contig chromosome;
    genomes from the main --genome_info file that are absent from it are submitted
    normally. The file needs columns genome_name, chromosome_name, chromosome_type
    and chromosome_topology; chromosome_location is optional. Values are
    not validated against the controlled vocabularies here - that happens per genome via
    validate_single_contig_fields().
    Args:
        single_contig_info (Path, optional): Path to the --single-contig-info TSV, or
            None/empty when the flag wasn't given.
    Returns:
        dict: genome_name -> {CHROMOSOME_NAME_FIELD: str, CHROMOSOME_TYPE_FIELD: str,
            CHROMOSOME_TOPOLOGY_FIELD: str, CHROMOSOME_LOCATION_FIELD: str}.
            Empty when single_contig_info was not given.
    Raises:
        ValueError: If a mandatory column is missing, or genome_name values are duplicated.
    """
    if not single_contig_info:
        return {}

    header_columns = pd.read_csv(single_contig_info, sep="\t", nrows=0).columns
    mandatory_columns = [
        GENOME_NAME_FIELD,
        CHROMOSOME_NAME_FIELD,
        CHROMOSOME_TYPE_FIELD,
        CHROMOSOME_TOPOLOGY_FIELD,
    ]
    missing_columns = [c for c in mandatory_columns if c not in header_columns]
    if missing_columns:
        raise ValueError(
            "The single-contig metadata file (--single-contig-info) is missing mandatory "
            f"column(s): {', '.join(missing_columns)}"
        )

    columns_to_read = list(mandatory_columns)
    if CHROMOSOME_LOCATION_FIELD in header_columns:
        columns_to_read.append(CHROMOSOME_LOCATION_FIELD)

    metadata = pd.read_csv(single_contig_info, sep="\t", usecols=columns_to_read)
    if CHROMOSOME_LOCATION_FIELD not in metadata.columns:
        metadata[CHROMOSOME_LOCATION_FIELD] = None

    if metadata[GENOME_NAME_FIELD].nunique() != metadata[GENOME_NAME_FIELD].size:
        raise ValueError("Duplicate genome_name values found in the single-contig metadata file (--single-contig-info)")

    return metadata.set_index(GENOME_NAME_FIELD).to_dict(orient="index")


def validate_single_contig_fields(genome_name: str, genome: dict) -> list:
    """
    Validates a single-contig genome's chromosome metadata.

    chromosome_name has no fixed vocabulary: it must be a digit string
    (chromosome/plasmid number) or "MIT" for the mitochondrial chromosome
    (CHROMOSOME_NAME_REGEX). chromosome_type and chromosome_topology
    must each be one of the corresponding constants vocabulary. chromosome_location
    is optional and only checked against its vocabulary when a value is present.

    The values are normalised in place (NaN/None -> "") so downstream writers
    receive clean strings.
    Args:
        genome_name (str): Name of the genome, used in error messages.
        genome (dict): The genome's metadata row; mutated in place.
    Returns:
        list: Human-readable error messages, empty when every field is valid.
    """
    name = normalise_na(genome.get(CHROMOSOME_NAME_FIELD))
    chromosome_type = normalise_na(genome.get(CHROMOSOME_TYPE_FIELD))
    topology = normalise_na(genome.get(CHROMOSOME_TOPOLOGY_FIELD))
    location = normalise_na(genome.get(CHROMOSOME_LOCATION_FIELD))

    genome[CHROMOSOME_NAME_FIELD] = name
    genome[CHROMOSOME_TYPE_FIELD] = chromosome_type
    genome[CHROMOSOME_TOPOLOGY_FIELD] = topology
    genome[CHROMOSOME_LOCATION_FIELD] = location

    errors = []
    if not CHROMOSOME_NAME_REGEX.match(name):
        errors.append(
            f"Genome '{genome_name}': {CHROMOSOME_NAME_FIELD} '{name}' must be a digit "
            "string (e.g. '1', '2', ...) or 'MIT' for the mitochondrial chromosome."
        )
    if chromosome_type not in CHROMOSOME_TYPES:
        errors.append(
            f"Genome '{genome_name}': {CHROMOSOME_TYPE_FIELD} '{chromosome_type}' is not one of "
            f"{sorted(CHROMOSOME_TYPES)}."
        )
    if topology not in CHROMOSOME_TOPOLOGIES:
        errors.append(
            f"Genome '{genome_name}': {CHROMOSOME_TOPOLOGY_FIELD} '{topology}' is not one of "
            f"{CHROMOSOME_TOPOLOGIES}."
        )
    if location and location not in CHROMOSOME_LOCATIONS_LIST:
        errors.append(
            f"Genome '{genome_name}': {CHROMOSOME_LOCATION_FIELD} '{location}' is not one of "
            f"{CHROMOSOME_LOCATIONS_LIST}."
        )
    return errors


def write_chromosome_list(
    manifest_dir: Path,
    alias: str,
    contig_id: str,
    chromosome_name: str,
    chromosome_type: str,
    chromosome_topology: str,
    chromosome_location: str,
) -> Path:
    """
    Writes a gzip-compressed chromosome list file for a single-contig genome.

    The file has one tab-separated line per replicon, here always a single line:
    ``<object_name>\t<chromosome_name>\t<topology>-<type>[\t<chromosome_location>]``
    e.g. ``contig_1\t1\tCircular-Chromosome`` or ``contig_2\tMIT\tLinear-Chromosome``
    (see https://ena-docs.readthedocs.io/en/latest/submit/fileprep/assembly.html#chromosome-list-file).

    Args:
        manifest_dir (Path): Directory the manifest files (and this chromosome list) live in.
        alias (str): Genome alias, used to name the output file.
        contig_id (str): Identifier of the genome's single contig (OBJECT_NAME,
            must match the FASTA header).
        chromosome_name (str): CHROMOSOME_NAME value: a digit string or "MIT".
        chromosome_type (str): chromosome type, one of constants.CHROMOSOME_TYPES.
        chromosome_topology (str): chromosome topology, one of
            constants.CHROMOSOME_TOPOLOGIES. Combined with chromosome_type
            (as "<topology>-<type>") to form the third column.
        chromosome_location (str): CHROMOSOME_LOCATION value; appended as a fourth
            column only when a non-empty value is given.
    Returns:
        Path: Path of the written chromosome list file.
    """
    # the chromosome values are validated and normalised upstream via
    # validate_single_contig_fields(); genomes with invalid values never reach this point
    chromosome_list_path = manifest_dir / f"{alias}_chromosome_list.txt.gz"
    line = f"{contig_id}\t{chromosome_name}\t{chromosome_topology}-{chromosome_type}"
    if chromosome_location:
        line += f"\t{chromosome_location}"
    with gzip.open(chromosome_list_path, "wt") as f:
        f.write(line + "\n")
    return chromosome_list_path


class SingleContigLog:
    """
    Dedicated log for problems found while processing --single-contig-info (invalid
    chromosome metadata, unmatched genome_name values, ...).

    The log file and its handler are created lazily, on the first write, so a run
    with no single-contig issues to report leaves no empty file behind.
    """

    def __init__(self, path: Path):
        """
        Args:
            path (Path): Where the log file will be written.
        """
        self.path = path
        self._logger = None

    def __str__(self) -> str:
        return str(self.path)

    def _get_logger(self) -> logging.Logger:
        if self._logger is None:
            sc_logger = logging.getLogger("genomeuploader.single_contig")
            sc_logger.setLevel(logging.INFO)
            sc_logger.propagate = False
            # drop stale handlers if this is instantiated more than once in a process
            # (e.g. one GenomeUpload per test), so each instance logs to its own path
            for handler in list(sc_logger.handlers):
                sc_logger.removeHandler(handler)
                handler.close()
            file_handler = logging.FileHandler(self.path, mode="w")
            file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
            sc_logger.addHandler(file_handler)
            self._logger = sc_logger
        return self._logger

    def error(self, message: str) -> None:
        """Records an error-level message."""
        self._get_logger().error(message)

    def warning(self, message: str) -> None:
        """Records a warning-level message."""
        self._get_logger().warning(message)

    def exists(self) -> bool:
        """Whether the log file has been created (i.e. something was logged)."""
        return self.path.exists()

    def read_text(self) -> str:
        """Reads back the log file's contents."""
        return self.path.read_text()
