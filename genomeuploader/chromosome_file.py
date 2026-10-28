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
Everything specific to submitting a genome as ENA chromosome-level genome assembly: parsing the
optional --chromosome-info TSV, validating its values and cross-checking them against
the genomes' fasta files up front, and writing the ENA chromosome list file.
See the ENA docs:
https://ena-docs.readthedocs.io/en/latest/submit/assembly/metagenome/mag.html#chromosome-assembly
"""

from dataclasses import dataclass
from fileinput import hook_compressed
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO

from genomeuploader.constants import (
    COMPRESSED_FASTA_SUFFIXES,
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

# columns every --chromosome-info file must have; chromosome_location is optional
MANDATORY_COLUMNS = [
    GENOME_NAME_FIELD,
    GENOME_CONTIG_FIELD,
    CHROMOSOME_NAME_FIELD,
    CHROMOSOME_TYPE_FIELD,
    CHROMOSOME_TOPOLOGY_FIELD,
]


@dataclass(frozen=True)
class ChromosomeRecord:
    """
    One row of --chromosome-info: a genome's contig submitted as a chromosome.

    Records are built with from_row(), which lowercases chromosome_type and
    chromosome_topology, and are immutable afterwards. Each field has its own
    validate_* method; validate() runs them all. contig_name is not checked here -
    validate_contig_names() needs the genome's fasta for that - and uniqueness of
    contig_name/chromosome_name within a genome is checked across rows by
    load_chromosome_metadata().
    """

    genome_name: str
    contig_name: str
    chromosome_name: str
    chromosome_type: str
    chromosome_topology: str
    chromosome_location: str = ""

    @classmethod
    def from_row(cls, row: dict) -> "ChromosomeRecord":
        """
        Builds a record from a --chromosome-info row as read by load_chromosome_metadata()
        (every value a stripped string, "" when empty).
        Args:
            row (dict): Column name -> cell value.
        Returns:
            ChromosomeRecord: The record, with chromosome_type and chromosome_topology lowercased.
        """
        return cls(
            genome_name=row[GENOME_NAME_FIELD],
            contig_name=row[GENOME_CONTIG_FIELD],
            chromosome_name=row[CHROMOSOME_NAME_FIELD],
            chromosome_type=row[CHROMOSOME_TYPE_FIELD].lower(),
            chromosome_topology=row[CHROMOSOME_TOPOLOGY_FIELD].lower(),
            chromosome_location=row[CHROMOSOME_LOCATION_FIELD],
        )

    def validate(self) -> list[str]:
        """
        Runs every field validation. Empty values are not reported here -
        load_chromosome_metadata() rejects empty mandatory fields - so each field is
        only checked when a value is present.
        Returns:
            list[str]: Human-readable error messages, empty when every field is valid.
        """
        return (
            self.validate_chromosome_name()
            + self.validate_chromosome_type()
            + self.validate_chromosome_topology()
            + self.validate_chromosome_location()
        )

    def validate_chromosome_name(self) -> list[str]:
        """
        chromosome_name has no fixed vocabulary; it appears as the /chromosome, /plasmid or
        /segment qualifier in the EMBL-Bank flat files. Names must:
        - match the pattern ^[A-Za-z0-9][A-Za-z0-9_#-.]*$, i.e. start with a letter or digit,
          then only letters, digits, '_', '#', '-' or '.' (CHROMOSOME_NAME_REGEX),
        - be shorter than 33 characters (CHROMOSOME_NAME_MAX_LENGTH),
        - be unique within an assembly (checked across rows by load_chromosome_metadata()),
        - not contain any of CHROMOSOME_NAME_FORBIDDEN_SUBSTRINGS as part of their name
          (case insensitive): 'chr', 'chrm', 'chrom', 'chromosome', 'linkage group',
          'linkage-group', 'linkage_group', 'plasmid'.
        Returns:
            list[str]: One message per broken rule, empty when the name is valid or empty.
        """
        name = self.chromosome_name
        if not name:
            return []
        errors = []
        if not CHROMOSOME_NAME_REGEX.match(name):
            errors.append(
                f"Genome '{self.genome_name}': {CHROMOSOME_NAME_FIELD} '{name}' must match the pattern "
                f"{CHROMOSOME_NAME_REGEX.pattern} (start with a letter or digit, then only letters, "
                "digits, '_', '#', '-' or '.')."
            )
        if len(name) > CHROMOSOME_NAME_MAX_LENGTH:
            errors.append(
                f"Genome '{self.genome_name}': {CHROMOSOME_NAME_FIELD} '{name}' must be shorter than "
                f"{CHROMOSOME_NAME_MAX_LENGTH + 1} characters (got {len(name)})."
            )
        forbidden = [s for s in CHROMOSOME_NAME_FORBIDDEN_SUBSTRINGS if s in name.lower()]
        if forbidden:
            errors.append(
                f"Genome '{self.genome_name}': {CHROMOSOME_NAME_FIELD} '{name}' must not contain "
                f"{', '.join(repr(s) for s in forbidden)} (case insensitive)."
            )
        return errors

    def validate_chromosome_type(self) -> list[str]:
        """
        chromosome_type must be one of CHROMOSOME_TYPES (matched case insensitively, as
        from_row() lowercases it).
        Returns:
            list[str]: An error message, or empty when the type is valid or empty.
        """
        if self.chromosome_type and self.chromosome_type not in CHROMOSOME_TYPES:
            return [
                f"Genome '{self.genome_name}': {CHROMOSOME_TYPE_FIELD} '{self.chromosome_type}' "
                f"is not one of {sorted(CHROMOSOME_TYPES)}."
            ]
        return []

    def validate_chromosome_topology(self) -> list[str]:
        """
        chromosome_topology must be one of CHROMOSOME_TOPOLOGIES (matched case
        insensitively, as from_row() lowercases it).
        Returns:
            list[str]: An error message, or empty when the topology is valid or empty.
        """
        if self.chromosome_topology and self.chromosome_topology not in CHROMOSOME_TOPOLOGIES:
            return [
                f"Genome '{self.genome_name}': {CHROMOSOME_TOPOLOGY_FIELD} '{self.chromosome_topology}' "
                f"is not one of {CHROMOSOME_TOPOLOGIES}."
            ]
        return []

    def validate_chromosome_location(self) -> list[str]:
        """
        chromosome_location is optional; when given it must be one of
        CHROMOSOME_LOCATIONS_LIST (exact match).
        Returns:
            list[str]: An error message, or empty when the location is valid or empty.
        """
        if self.chromosome_location and self.chromosome_location not in CHROMOSOME_LOCATIONS_LIST:
            return [
                f"Genome '{self.genome_name}': {CHROMOSOME_LOCATION_FIELD} '{self.chromosome_location}' "
                f"is not one of {CHROMOSOME_LOCATIONS_LIST}."
            ]
        return []

    def to_chromosome_list_line(self) -> str:
        """
        Formats the record as a chromosome list file line (tab-separated):
        ``<contig_name>\t<chromosome_name>\t<topology>-<type>[\t<chromosome_location>]``,
        with chromosome_location only appended when non-empty.
        Returns:
            str: The line, without a trailing newline.
        """
        line = f"{self.contig_name}\t{self.chromosome_name}\t{self.chromosome_topology}-{self.chromosome_type}"
        if self.chromosome_location:
            line += f"\t{self.chromosome_location}"
        return line


def get_fasta_contig_ids(fasta_path: Path) -> list[str]:
    """
    Extracts sequence identifiers from a gzip- or bzip2-compressed FASTA file.
    Args:
        fasta_path (Path): Path to the compressed FASTA file (.gz or .bz2).
    Returns:
        list[str]: Contig/sequence identifiers found in the FASTA file.
    Raises:
        ValueError: If the file doesn't have a .gz/.bz2 extension or its content isn't
            compressed in the format its extension says.
    """
    # hook_compressed picks the decompressor from the extension and would read any
    # other file as plain text, so the extension is checked first
    if Path(fasta_path).suffix not in COMPRESSED_FASTA_SUFFIXES:
        raise ValueError(
            f"Genome file {fasta_path} must be compressed with gzip or bzip2 "
            f"(extension one of {COMPRESSED_FASTA_SUFFIXES})."
        )
    with hook_compressed(fasta_path, "r") as contigs:
        try:
            return [record.id for record in SeqIO.parse(contigs, "fasta")]
        # gzip.BadGzipFile and bz2's "Invalid data stream" are both OSErrors
        except OSError:
            raise ValueError(
                f"Genome file {fasta_path} is not valid {Path(fasta_path).suffix} compressed data."
            )


def read_chromosome_info(chromosome_info: Path) -> list[dict]:
    """
    Reads the --chromosome-info TSV into cleaned rows, without validating their values.

    Checks that every column in MANDATORY_COLUMNS is present, then reads only those
    columns plus chromosome_location when the file has it. Every cell is read as a
    string with surrounding whitespace stripped, and empty cells are "" (never NaN).
    A missing chromosome_location column is added with "" in every row, so each
    returned row has the same keys.
    Args:
        chromosome_info (Path): Path to the --chromosome-info TSV.
    Returns:
        list[dict]: One dict per data row, mapping column name -> stripped string value,
            in file order (the first data row is line 2 of the file).
    Raises:
        ValueError: If any mandatory column is missing from the header.
    """
    header_columns = pd.read_csv(chromosome_info, sep="\t", nrows=0).columns
    missing_columns = [c for c in MANDATORY_COLUMNS if c not in header_columns]
    if missing_columns:
        raise ValueError(
            "The chromosome metadata file (--chromosome-info) is missing mandatory "
            f"column(s): {', '.join(missing_columns)}"
        )

    columns_to_read = list(MANDATORY_COLUMNS)
    if CHROMOSOME_LOCATION_FIELD in header_columns:
        columns_to_read.append(CHROMOSOME_LOCATION_FIELD)

    # read everything as str with empty cells as "" (not NaN): otherwise a column with an
    # empty cell is parsed as float and e.g. chromosome_name "1" becomes "1.0"
    metadata = pd.read_csv(chromosome_info, sep="\t", dtype=str, keep_default_na=False, usecols=columns_to_read)
    metadata = metadata.apply(lambda col: col.str.strip())
    if CHROMOSOME_LOCATION_FIELD not in metadata.columns:
        metadata[CHROMOSOME_LOCATION_FIELD] = ""

    return metadata.to_dict(orient="records")


def validate_chromosome_records(rows: list[dict], records: list[ChromosomeRecord]) -> list[str]:
    """
    Validates every --chromosome-info row and collects all problems instead of stopping
    at the first one.

    Checks that:
    - every mandatory field (MANDATORY_COLUMNS) is non-empty in every row, reported with
      the row's line number in the file so it can be found even when genome_name is empty,
    - every record's values are valid (see ChromosomeRecord.validate()),
    - contig_name and chromosome_name are each unique within a genome_name, as each
      genome is one assembly (empty values are skipped, as they're already reported).
    Whether each contig_name exists in the genome's fasta is not checked here - see
    validate_contig_names().
    Args:
        rows (list[dict]): Cleaned rows, as returned by read_chromosome_info().
        records (list[ChromosomeRecord]): The records built from those rows, in the same order.
    Returns:
        list[str]: Human-readable error messages, empty when every row is valid.
    """
    validation_errors = []

    # row numbers are 1-based file lines (line 1 is the header), so they still identify
    # the row when genome_name itself is the empty field
    for row_number, (row, record) in enumerate(zip(rows, records), start=2):
        empty_fields = [c for c in MANDATORY_COLUMNS if not row[c]]
        if empty_fields:
            validation_errors.append(
                f"Row {row_number} (genome '{record.genome_name}'): mandatory field(s) "
                f"{', '.join(empty_fields)} must not be empty."
            )
        validation_errors.extend(record.validate())

    # different chromosome records for the same genome (one assembly) must point at
    # different contigs and have different chromosome names
    for field, requirement in (
        (GENOME_CONTIG_FIELD, "each record must point at a different contig"),
        (CHROMOSOME_NAME_FIELD, "chromosome names must be unique within an assembly"),
    ):
        seen_pairs = set()
        duplicate_pairs = set()
        for record in records:
            # the column name constants are also the ChromosomeRecord attribute names
            pair = (record.genome_name, getattr(record, field))
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
    return validation_errors


def load_chromosome_metadata(chromosome_info: Optional[Path]) -> dict[str, list[ChromosomeRecord]]:
    """
    Loads and validates the optional chromosome metadata TSV (--chromosome-info) and
    groups its rows into ChromosomeRecords per genome.

    A genome can have one or more chromosome records - one row per contig submitted as
    a chromosome (e.g. chr1, chr2, a plasmid, the mitochondrial contig, ...), all
    sharing the same genome_name. The file is read by read_chromosome_info() (column
    presence, whitespace stripping), each row is turned into a ChromosomeRecord with
    ChromosomeRecord.from_row(), and all rows are checked by validate_chromosome_records();
    every problem found is raised together in one error. Whether each contig_name exists
    in the genome's fasta is checked separately by validate_contig_names(), once genome
    paths are known. This either returns fully valid data or raises - callers never need
    to handle invalid rows.
    Args:
        chromosome_info (Path, optional): Path to the --chromosome-info TSV, or
            None/empty when the flag wasn't given.
    Returns:
        dict[str, list[ChromosomeRecord]]: genome_name -> that genome's chromosome
            records, one per row. Empty when chromosome_info was not given.
    Raises:
        ValueError: If a mandatory column is missing (raised by read_chromosome_info()),
            or validate_chromosome_records() reports any problem.
    """
    if not chromosome_info:
        return {}

    rows = read_chromosome_info(chromosome_info)
    records = [ChromosomeRecord.from_row(row) for row in rows]

    validation_errors = validate_chromosome_records(rows, records)

    if validation_errors:
        raise ValueError(
            "The chromosome metadata file (--chromosome-info) failed validation:\n- "
            + "\n- ".join(validation_errors)
        )

    chromosome_metadata = {}
    for record in records:
        chromosome_metadata.setdefault(record.genome_name, []).append(record)

    return chromosome_metadata


def validate_contig_names(chromosome_metadata: dict[str, list[ChromosomeRecord]], genome_paths: dict) -> None:
    """
    Checks that every contig_name in chromosome_metadata exists in its genome's fasta file.

    Only genome_name values present in genome_paths (i.e. genomes in the current
    --genome_info batch) are checked; a genome_name from --chromosome-info that doesn't
    match any genome is reported separately by the caller (extract_genomes_info()).
    Args:
        chromosome_metadata (dict[str, list[ChromosomeRecord]]): genome_name -> chromosome
            records, as returned by load_chromosome_metadata().
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
            contig_name = record.contig_name
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


def write_chromosome_list(manifest_dir: Path, alias: str, chromosome_records: list[ChromosomeRecord]) -> Path:
    """
    Writes a chromosome list file for a genome, one line per chromosome record.

    Each line is formatted by ChromosomeRecord.to_chromosome_list_line(),
    e.g. ``chr1\t1\tlinear-chromosome`` or ``chrMT\tMIT\tcircular-chromosome\tMitochondrion``
    (see https://ena-docs.readthedocs.io/en/latest/submit/fileprep/assembly.html#chromosome-list-file).

    Args:
        manifest_dir (Path): Directory the manifest files (and this chromosome list) live in.
        alias (str): Genome alias, used to name the output file.
        chromosome_records (list[ChromosomeRecord]): The genome's chromosome records
            (genome_info[gen]["chromosome_records"]).
    Returns:
        Path: Path of the written chromosome list file.
    """
    # the chromosome values are validated and normalised upstream via
    # load_chromosome_metadata(); genomes with invalid values never reach this point
    chromosome_list_path = manifest_dir / f"{alias}_chromosome_list.txt"
    lines = [record.to_chromosome_list_line() for record in chromosome_records]
    with open(chromosome_list_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return chromosome_list_path
