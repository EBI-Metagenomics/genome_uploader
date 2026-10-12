#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2017-2026 EMBL - European Bioinformatics Institute
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
Everything specific to submitting a genome's contigs as ENA chromosome-level assembly: parsing the
optional --chromosome-info TSV, validating its values and cross-checking them against
the genomes' fasta files up front, and writing the ENA chromosome list file.
See the ENA docs:
https://ena-docs.readthedocs.io/en/latest/submit/assembly/metagenome/mag.html#chromosome-assembly
"""

import gzip
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO

from genomeuploader.constants import (
    GENOME_NAME_FIELD,
    GENOME_CONTIG_FIELD,
    CHROMOSOME_NAME_FIELD,
    CHROMOSOME_NAME_REGEX,
    CHROMOSOME_NAME_MAX_LENGTH,
    CHROMOSOME_NAME_FORBIDDEN_SUBSTRINGS,
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


def get_fasta_contig_ids(fasta_path: Path) -> list[str]:
    """
    Extracts sequence identifiers from a gzip-compressed FASTA file.
    Args:
        fasta_path (Path): Path to the gzip-compressed FASTA file.
    Returns:
        list[str]: Contig/sequence identifiers found in the FASTA file.
    """
    with gzip.open(fasta_path, "rt") as contigs:
        try:
            return [record.id for record in SeqIO.parse(contigs, "fasta")]
        except gzip.BadGzipFile:
            raise ValueError(f"Genome file {fasta_path} is not gzip-compressed. Genome files must be gzip-compressed.")


def load_chromosome_metadata(chromosome_info: Optional[Path]) -> dict:
    """
    Loads, validates and normalises the optional chromosome metadata TSV (--chromosome-info).

    A genome can have one or more chromosome records - one row per contig submitted as
    a chromosome (e.g. chr1, chr2, a plasmid, the mitochondrial contig, ...), all
    sharing the same genome_name. Each row needs columns genome_name, contig_name
    (identifies which of that genome's fasta sequences the row describes - existence
    is checked separately by validate_contig_names(), once genome paths are known),
    chromosome_name, chromosome_type and chromosome_topology; chromosome_location is
    optional. Every mandatory field must be non-empty in every row, every row's
    chromosome metadata is validated against the controlled vocabularies
    (see validate_chromosome_fields()) as part of loading, and both contig_name and
    chromosome_name must be unique within each genome_name (two chromosome records for
    the same genome can't point at the same contig or share a chromosome name, as each
    genome is one assembly), so this either returns fully valid, normalised
    data or raises - callers never need to handle invalid rows.
    Args:
        chromosome_info (Path, optional): Path to the --chromosome-info TSV, or
            None/empty when the flag wasn't given.
    Returns:
        dict: genome_name -> list of dicts {GENOME_CONTIG_FIELD: str, CHROMOSOME_NAME_FIELD: str,
            CHROMOSOME_TYPE_FIELD: str, CHROMOSOME_TOPOLOGY_FIELD: str,
            CHROMOSOME_LOCATION_FIELD: str}, one entry per chromosome record for that
            genome. Empty when chromosome_info was not given.
    Raises:
        ValueError: If a mandatory column is missing, a mandatory field is empty in any
            row, the same genome_name/contig_name or genome_name/chromosome_name pair
            appears more than once, or any row's chromosome metadata fails validation.
    """
    if not chromosome_info:
        return {}

    header_columns = pd.read_csv(chromosome_info, sep="\t", nrows=0).columns
    mandatory_columns = [
        GENOME_NAME_FIELD,
        GENOME_CONTIG_FIELD,
        CHROMOSOME_NAME_FIELD,
        CHROMOSOME_TYPE_FIELD,
        CHROMOSOME_TOPOLOGY_FIELD,
    ]
    missing_columns = [c for c in mandatory_columns if c not in header_columns]
    if missing_columns:
        raise ValueError(
            "The chromosome metadata file (--chromosome-info) is missing mandatory "
            f"column(s): {', '.join(missing_columns)}"
        )

    columns_to_read = list(mandatory_columns)
    if CHROMOSOME_LOCATION_FIELD in header_columns:
        columns_to_read.append(CHROMOSOME_LOCATION_FIELD)

    # read everything as str: otherwise a column with an empty cell is parsed as float
    # and e.g. chromosome_name "1" becomes "1.0"
    metadata = pd.read_csv(chromosome_info, sep="\t", usecols=columns_to_read, dtype=str)
    if CHROMOSOME_LOCATION_FIELD not in metadata.columns:
        metadata[CHROMOSOME_LOCATION_FIELD] = None

    records = metadata.to_dict(orient="records")
    for record in records:
        record[GENOME_NAME_FIELD] = normalise_na(record.get(GENOME_NAME_FIELD))
        record[GENOME_CONTIG_FIELD] = normalise_na(record.get(GENOME_CONTIG_FIELD))

    validation_errors = []

    # row numbers are 1-based file lines (line 1 is the header), so they still identify
    # the row when genome_name itself is the empty field
    for row_number, record in enumerate(records, start=2):
        empty_fields = [c for c in mandatory_columns if not normalise_na(record.get(c))]
        if empty_fields:
            validation_errors.append(
                f"Row {row_number} (genome '{record[GENOME_NAME_FIELD]}'): mandatory field(s) "
                f"{', '.join(empty_fields)} must not be empty."
            )
        validation_errors.extend(validate_chromosome_fields(record[GENOME_NAME_FIELD], record))

    # different chromosome records for the same genome (one assembly) must point at
    # different contigs and have different chromosome names
    for field, requirement in (
        (GENOME_CONTIG_FIELD, "each record must point at a different contig"),
        (CHROMOSOME_NAME_FIELD, "chromosome names must be unique within an assembly"),
    ):
        seen_pairs = set()
        duplicate_pairs = set()
        for record in records:
            pair = (record[GENOME_NAME_FIELD], record[field])
            # empty values are already reported as missing mandatory fields
            if not pair[1]:
                continue
            if pair in seen_pairs:
                duplicate_pairs.add(pair)
            seen_pairs.add(pair)
        for genome_name, value in sorted(duplicate_pairs):
            validation_errors.append(
                f"Genome '{genome_name}': {field} '{value}' is listed in more than "
                f"one chromosome record; {requirement}."
            )

    if validation_errors:
        raise ValueError(
            "The chromosome metadata file (--chromosome-info) failed validation:\n- "
            + "\n- ".join(validation_errors)
        )

    chromosome_metadata = {}
    for record in records:
        chromosome_metadata.setdefault(record[GENOME_NAME_FIELD], []).append(record)

    return chromosome_metadata


def validate_contig_names(chromosome_metadata: dict, genome_paths: dict) -> None:
    """
    Checks that every contig_name in chromosome_metadata exists in its genome's fasta file.

    Only genome_name values present in genome_paths (i.e. genomes in the current
    --genome_info batch) are checked; a genome_name from --chromosome-info that doesn't
    match any genome is a separate, non-fatal concern handled by the caller.
    Args:
        chromosome_metadata (dict): genome_name -> list of chromosome records, as
            returned by load_chromosome_metadata().
        genome_paths (dict): genome_name -> Path to that genome's fasta file.
    Raises:
        ValueError: If any contig_name is not found among its genome's fasta sequence
            headers.
    """
    errors = []
    for genome_name, records in chromosome_metadata.items():
        genome_path = genome_paths.get(genome_name)
        if genome_path is None:
            continue
        contig_ids = set(get_fasta_contig_ids(genome_path))
        for record in records:
            contig_name = record[GENOME_CONTIG_FIELD]
            if contig_name not in contig_ids:
                errors.append(
                    f"Genome '{genome_name}': {GENOME_CONTIG_FIELD} '{contig_name}' was not found "
                    f"among the sequence headers in {genome_path}."
                )
    if errors:
        raise ValueError(
            "The chromosome metadata file (--chromosome-info) references contig(s) that don't "
            "exist in the corresponding genome fasta file(s):\n- " + "\n- ".join(errors)
        )


def validate_chromosome_fields(genome_name: str, genome: dict) -> list:
    """
    Validates a single chromosome record's metadata (everything but contig_name, which
    is checked separately by validate_contig_names() since it needs the genome's fasta).

    chromosome_name has no fixed vocabulary; it appears as the /chromosome, /plasmid or
    /segment qualifier in the EMBL-Bank flat files. Names must:
    - match the pattern ^[A-Za-z0-9][A-Za-z0-9_#-.]*$, i.e. start with a letter or digit,
      then only letters, digits, '_', '#', '-' or '.' (CHROMOSOME_NAME_REGEX),
    - be shorter than 33 characters (CHROMOSOME_NAME_MAX_LENGTH),
    - be unique within an assembly (checked across rows by load_chromosome_metadata()),
    - not contain any of CHROMOSOME_NAME_FORBIDDEN_SUBSTRINGS as part of their name
      (case insensitive): 'chr', 'chrm', 'chrom', 'chromosome', 'linkage group',
      'linkage-group', 'linkage_group', 'plasmid'.
    chromosome_type and chromosome_topology
    must each be one of the corresponding constants vocabulary; they are matched case
    insensitively and stored lowercased (e.g. "Circular" -> "circular"). chromosome_location
    is optional. Empty values are not reported here - load_chromosome_metadata()
    rejects empty mandatory fields - so every field is only checked against its
    vocabulary when a value is present.

    The values are normalised in place (NaN/None -> "") so downstream writers
    receive clean strings.
    Args:
        genome_name (str): Name of the genome, used in error messages.
        genome (dict): The chromosome record; mutated in place.
    Returns:
        list: Human-readable error messages, empty when every field is valid.
    """
    name = normalise_na(genome.get(CHROMOSOME_NAME_FIELD))
    chromosome_type = normalise_na(genome.get(CHROMOSOME_TYPE_FIELD)).lower()
    topology = normalise_na(genome.get(CHROMOSOME_TOPOLOGY_FIELD)).lower()
    location = normalise_na(genome.get(CHROMOSOME_LOCATION_FIELD))

    genome[CHROMOSOME_NAME_FIELD] = name
    genome[CHROMOSOME_TYPE_FIELD] = chromosome_type
    genome[CHROMOSOME_TOPOLOGY_FIELD] = topology
    genome[CHROMOSOME_LOCATION_FIELD] = location

    # empty mandatory values are reported by load_chromosome_metadata(), so only
    # non-empty values are checked against their vocabularies here
    errors = []
    if name:
        if not CHROMOSOME_NAME_REGEX.match(name):
            errors.append(
                f"Genome '{genome_name}': {CHROMOSOME_NAME_FIELD} '{name}' must match the pattern "
                f"{CHROMOSOME_NAME_REGEX.pattern} (start with a letter or digit, then only letters, "
                "digits, '_', '#', '-' or '.')."
            )
        if len(name) > CHROMOSOME_NAME_MAX_LENGTH:
            errors.append(
                f"Genome '{genome_name}': {CHROMOSOME_NAME_FIELD} '{name}' must be shorter than "
                f"{CHROMOSOME_NAME_MAX_LENGTH + 1} characters (got {len(name)})."
            )
        forbidden = [s for s in CHROMOSOME_NAME_FORBIDDEN_SUBSTRINGS if s in name.lower()]
        if forbidden:
            errors.append(
                f"Genome '{genome_name}': {CHROMOSOME_NAME_FIELD} '{name}' must not contain "
                f"{', '.join(repr(s) for s in forbidden)} (case insensitive)."
            )
    if chromosome_type and chromosome_type not in CHROMOSOME_TYPES:
        errors.append(
            f"Genome '{genome_name}': {CHROMOSOME_TYPE_FIELD} '{chromosome_type}' is not one of "
            f"{sorted(CHROMOSOME_TYPES)}."
        )
    if topology and topology not in CHROMOSOME_TOPOLOGIES:
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


def write_chromosome_list(manifest_dir: Path, alias: str, chromosome_records: list) -> Path:
    """
    Writes a chromosome list file for a genome, one line per chromosome record.

    Each line is tab-separated:
    ``<contig_name>\t<chromosome_name>\t<topology>-<type>[\t<chromosome_location>]``
    e.g. ``chr1\t1\tcircular-chromosome`` or ``chrMT\tMIT\tlinear-chromosome``
    (see https://ena-docs.readthedocs.io/en/latest/submit/fileprep/assembly.html#chromosome-list-file).

    Args:
        manifest_dir (Path): Directory the manifest files (and this chromosome list) live in.
        alias (str): Genome alias, used to name the output file.
        chromosome_records (list): The genome's chromosome records (genome_info[gen]
            ["chromosome_records"]), each a dict with GENOME_CONTIG_FIELD (OBJECT_NAME,
            must match a FASTA header), CHROMOSOME_NAME_FIELD (e.g. "1" or "MIT"),
            CHROMOSOME_TYPE_FIELD, CHROMOSOME_TOPOLOGY_FIELD, and CHROMOSOME_LOCATION_FIELD
            (appended as a fourth column only when non-empty).
    Returns:
        Path: Path of the written chromosome list file.
    """
    # the chromosome values are validated and normalised upstream via
    # load_chromosome_metadata(); genomes with invalid values never reach this point
    chromosome_list_path = manifest_dir / f"{alias}_chromosome_list.txt"
    lines = []
    for record in chromosome_records:
        line = (
            f"{record[GENOME_CONTIG_FIELD]}\t{record[CHROMOSOME_NAME_FIELD]}\t"
            f"{record[CHROMOSOME_TOPOLOGY_FIELD]}-{record[CHROMOSOME_TYPE_FIELD]}"
        )
        if record[CHROMOSOME_LOCATION_FIELD]:
            line += f"\t{record[CHROMOSOME_LOCATION_FIELD]}"
        lines.append(line)
    with open(chromosome_list_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return chromosome_list_path
