import gzip
from pathlib import Path

import pytest

from genomeuploader.genome_upload import *  # noqa: F401,F403 - GenomeUpload, compute_manifests
from genomeuploader.constants import (
    GENOME_CONTIG_FIELD,
    CHROMOSOME_NAME_FIELD,
    CHROMOSOME_TYPE_FIELD,
    CHROMOSOME_TOPOLOGY_FIELD,
    CHROMOSOME_LOCATION_FIELD,
)

# All input files (genome metadata tsv, chromosome-info tsv, fasta.gz) used by these
# tests are built on the fly under `tmp_path` - nothing here is read from tests/fixtures.

GENOME_TSV_HEADER = [
    "genome_name",
    "genome_path",
    "accessions",
    "assembly_software",
    "binning_software",
    "binning_parameters",
    "stats_generation_software",
    "completeness",
    "contamination",
    "genome_coverage",
    "metagenome",
    "co-assembly",
    "broad_environment",
    "local_environment",
    "environmental_medium",
    "rRNA_presence",
    "NCBI_lineage",
]

NCBI_LINEAGE = (
    "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;"
    "g__Lactobacillus;s__Lactobacillus crispatus"
)

CHROMOSOME_INFO_HEADER = [
    "genome_name",
    "contig_name",
    "chromosome_name",
    "chromosome_type",
    "chromosome_topology",
]


def _write_tsv(path: Path, header: list, rows: list) -> Path:
    lines = ["\t".join(header)] + ["\t".join(str(value) for value in row) for row in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def _write_fasta_gz(path: Path, contig_ids: list) -> Path:
    with gzip.open(path, "wt") as f:
        for contig_id in contig_ids:
            f.write(f">{contig_id}\nACGTACGTACGTACGTACGTA\n")
    return path


def _genome_row(genome_name: str, genome_path: Path, accessions: str = "ERR6769700") -> list:
    return [
        genome_name,
        str(genome_path),
        accessions,
        "megahit_v1.2.9",
        "MGnify-genomes-generation-pipeline_v1.0.0",
        "default",
        "CheckM2_v1.0.1",
        "90.81314",
        "0.59",
        "14.2",
        "chicken gut metagenome",
        "False",
        "chicken",
        "gut",
        "mucosa",
        "True",
        NCBI_LINEAGE,
    ]


def _write_genome_tsv(path: Path, rows: list) -> Path:
    return _write_tsv(path, GENOME_TSV_HEADER, rows)


def _write_chromosome_info_tsv(path: Path, rows: list) -> Path:
    return _write_tsv(path, CHROMOSOME_INFO_HEADER, rows)


def _base_args(tmp_path: Path, genome_info: Path, chromosome_info: Path = None, test_suffix: str = "chromosome-unittest") -> dict:
    args = {
        "bins": True,
        "live": False,
        "private": False,
        "tpa": False,
        "centre_name": "EMG",
        "force": False,
        "out": str(tmp_path),
        "upload_study": "ERP000000",
        "genome_info": str(genome_info),
        "test_suffix": test_suffix,
    }
    if chromosome_info is not None:
        args["chromosome_info"] = str(chromosome_info)
    return args


class Tests:
    def test_chromosome_genome_writes_chromosome_list(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "chromosome_bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path)])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["ERR6769700_bin.1", "contig1", "1", "Chromosome", "Circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path)
        gu = GenomeUpload(args)
        # the chromosome-info file is loaded and validated in __init__ already
        assert gu.chromosome_metadata["ERR6769700_bin.1"][0][CHROMOSOME_NAME_FIELD] == "1"

        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        assert genome_info[alias]["has_chromosome"] is True
        records = genome_info[alias]["chromosome_records"]
        assert len(records) == 1
        assert records[0][GENOME_CONTIG_FIELD] == "contig1"
        assert records[0][CHROMOSOME_NAME_FIELD] == "1"
        assert records[0][CHROMOSOME_TYPE_FIELD] == "Chromosome"
        assert records[0][CHROMOSOME_TOPOLOGY_FIELD] == "Circular"
        assert records[0][CHROMOSOME_LOCATION_FIELD] == ""

        # don't call API again, force assign info for the next step
        genome_info[alias]["run_ref"] = "ERR6769700"
        genome_info[alias]["study"] = "ERP000000"
        genome_info[alias]["sequencingMethod"] = "ILLUMINA"

        manifest_info = compute_manifests(genome_info)[alias]
        gu.generate_genome_manifest(manifest_info, {alias: "ERS0000001"})

        manifest_text = (gu.manifest_dir / f"{alias}.manifest").read_text()
        chromosome_list_path = gu.manifest_dir / f"{alias}_chromosome_list.txt.gz"

        assert chromosome_list_path.exists()
        assert f"CHROMOSOME_LIST\t{chromosome_list_path.resolve()}" in manifest_text
        with gzip.open(chromosome_list_path, "rt") as f:
            assert f.read() == "contig1\t1\tCircular-Chromosome\n"

    def test_genome_with_multiple_chromosome_records_writes_all_lines(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "chromosome_bin.fa.gz", ["chr1", "chr2"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path)])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [
                ["ERR6769700_bin.1", "chr1", "1", "Chromosome", "Circular"],
                ["ERR6769700_bin.1", "chr2", "2", "Plasmid", "Linear"],
            ],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path, test_suffix="chromosome-multi-unittest")
        gu = GenomeUpload(args)
        assert len(gu.chromosome_metadata["ERR6769700_bin.1"]) == 2

        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        assert genome_info[alias]["has_chromosome"] is True
        assert len(genome_info[alias]["chromosome_records"]) == 2

        genome_info[alias]["run_ref"] = "ERR6769700"
        genome_info[alias]["study"] = "ERP000000"
        genome_info[alias]["sequencingMethod"] = "ILLUMINA"

        manifest_info = compute_manifests(genome_info)[alias]
        gu.generate_genome_manifest(manifest_info, {alias: "ERS0000001"})

        chromosome_list_path = gu.manifest_dir / f"{alias}_chromosome_list.txt.gz"
        with gzip.open(chromosome_list_path, "rt") as f:
            assert f.read() == "chr1\t1\tCircular-Chromosome\nchr2\t2\tLinear-Plasmid\n"

    def test_no_chromosome_info_provided_defaults_to_false(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("MAG1", genome_path)])
        args = _base_args(tmp_path, genome_info_path, test_suffix="no-chromosome-info-unittest")
        gu = GenomeUpload(args)
        assert gu.chromosome_metadata == {}
        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        assert genome_info[alias]["has_chromosome"] is False
        assert genome_info[alias]["chromosome_records"] == []

    def test_genome_not_listed_in_chromosome_info_defaults_to_false(self, tmp_path, caplog):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("MAG1", genome_path)])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["genome_that_does_not_exist", "contig1", "1", "Chromosome", "Circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path, test_suffix="chromosome-unmatched-unittest")
        gu = GenomeUpload(args)

        with caplog.at_level("WARNING", logger="genomeuploader.genome_upload"):
            genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        # the only genome isn't listed in the chromosome info file, so it's submitted normally,
        # and the unmatched genome_name's contig_name is never checked against any fasta
        assert genome_info[alias]["has_chromosome"] is False

        assert "genome_that_does_not_exist" in caplog.text
        assert "do not match any genome" in caplog.text

    def test_chromosome_info_contig_name_not_in_fasta_raises(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "chromosome_bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path)])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["ERR6769700_bin.1", "contig_that_does_not_exist", "1", "Chromosome", "Circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path, test_suffix="chromosome-missing-contig-unittest")
        gu = GenomeUpload(args)

        # structural/vocabulary validation passes in __init__ - the fasta cross-check only
        # happens once genome paths are known, in extract_genomes_info()
        with pytest.raises(ValueError, match="was not found among the sequence headers"):
            gu.extract_genomes_info()

    def test_chromosome_info_duplicate_contig_for_genome_raises(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [
                ["ERR6769700_bin.1", "contig1", "1", "Chromosome", "Circular"],
                ["ERR6769700_bin.1", "contig1", "2", "Plasmid", "Linear"],
            ],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-duplicate-contig-unittest")

        with pytest.raises(ValueError, match="is listed in more than one chromosome record"):
            GenomeUpload(args)

    def test_chromosome_info_with_wrong_type_raises(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [
                ["ERR6769700_bin.1", "contig1", "1", "not_a_real_type", "Circular"],
                ["ERR6769700_bin.2", "contig1", "2", "Segment", "Circular"],
            ],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-wrong-type-unittest")

        # both rows carry a chromosome_type that is not one of CHROMOSOME_TYPES, so the
        # whole file is rejected before anything else (--genome_info doesn't even need to exist)
        with pytest.raises(ValueError, match="failed validation") as exc_info:
            GenomeUpload(args)

        assert "chromosome_type 'not_a_real_type' is not one of" in str(exc_info.value)
        assert "chromosome_type 'Segment' is not one of" in str(exc_info.value)
        assert "ERR6769700_bin.1" in str(exc_info.value)
        assert "ERR6769700_bin.2" in str(exc_info.value)

    def test_chromosome_info_with_wrong_name_format_raises(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["ERR6769700_bin.1", "contig1", "chr1", "Chromosome", "Circular"]],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-wrong-name-unittest")

        # chromosome_name "chr1" is neither a digit string nor "MIT"
        with pytest.raises(ValueError, match="must be a digit string"):
            GenomeUpload(args)

    def test_chromosome_info_empty_contig_name_raises(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["ERR6769700_bin.1", "", "1", "Chromosome", "Circular"]],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-empty-contig-unittest")

        with pytest.raises(ValueError, match="contig_name must not be empty"):
            GenomeUpload(args)

    def test_mitochondrion_name_writes_chromosome_list(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "chromosome_bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path)])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["ERR6769700_bin.1", "contig1", "MIT", "Chromosome", "Linear"]],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path, test_suffix="chromosome-mitochondrion-unittest")
        gu = GenomeUpload(args)
        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        records = genome_info[alias]["chromosome_records"]
        assert genome_info[alias]["has_chromosome"] is True
        assert records[0][CHROMOSOME_NAME_FIELD] == "MIT"
        assert records[0][CHROMOSOME_TYPE_FIELD] == "Chromosome"
        assert records[0][CHROMOSOME_TOPOLOGY_FIELD] == "Linear"

        genome_info[alias]["run_ref"] = "ERR6769700"
        genome_info[alias]["study"] = "ERP000000"
        genome_info[alias]["sequencingMethod"] = "ILLUMINA"

        manifest_info = compute_manifests(genome_info)[alias]
        gu.generate_genome_manifest(manifest_info, {alias: "ERS0000001"})

        chromosome_list_path = gu.manifest_dir / f"{alias}_chromosome_list.txt.gz"
        with gzip.open(chromosome_list_path, "rt") as f:
            assert f.read() == "contig1\tMIT\tLinear-Chromosome\n"

    def test_chromosome_info_missing_mandatory_column_raises(self, tmp_path):
        # chromosome_topology is missing, even though the others are present
        chromosome_info_path = tmp_path / "chromosome_info.tsv"
        chromosome_info_path.write_text(
            "genome_name\tcontig_name\tchromosome_name\tchromosome_type\nMAG1\tcontig1\t1\tChromosome\n"
        )

        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-missing-column-unittest")

        with pytest.raises(ValueError, match="chromosome_topology"):
            GenomeUpload(args)

    def test_chromosome_info_missing_contig_name_column_raises(self, tmp_path):
        chromosome_info_path = tmp_path / "chromosome_info.tsv"
        chromosome_info_path.write_text(
            "genome_name\tchromosome_name\tchromosome_type\tchromosome_topology\nMAG1\t1\tChromosome\tCircular\n"
        )

        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-missing-contig-column-unittest")

        with pytest.raises(ValueError, match="contig_name"):
            GenomeUpload(args)
