import gzip
from pathlib import Path

import pytest

from genomeuploader.genome_upload import *  # noqa: F401,F403 - GenomeUpload, compute_manifests, SINGLE_CONTIG_* constants

# All input files (genome metadata tsv, single-contig-info tsv, fasta.gz) used by these
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

SINGLE_CONTIG_INFO_HEADER = [
    "genome_name",
    "single_contig_name",
    "single_contig_type",
    "single_contig_topology",
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


def _write_single_contig_info_tsv(path: Path, rows: list) -> Path:
    return _write_tsv(path, SINGLE_CONTIG_INFO_HEADER, rows)


def _base_args(tmp_path: Path, genome_info: Path, single_contig_info: Path = None, test_suffix: str = "single-contig-unittest") -> dict:
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
    if single_contig_info is not None:
        args["single_contig_info"] = str(single_contig_info)
    return args


class Tests:
    def test_single_contig_genome_extracts_contig_id_and_writes_chromosome_list(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "single_contig_bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path)])
        single_contig_info_path = _write_single_contig_info_tsv(
            tmp_path / "single_contig_info.tsv",
            [["ERR6769700_bin.1", "1", "Chromosome", "Circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, single_contig_info_path)
        gu = GenomeUpload(args)
        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        assert genome_info[alias]["single_contig"] is True
        assert genome_info[alias]["contig_id"] == "contig1"
        # the chromosome columns are normalised in place and kept on the genome dict
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_NAME] == "1"
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_TYPE] == "Chromosome"
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_TOPOLOGY] == "Circular"
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_LOCATION] == ""
        # no chromosome metadata problems -> dedicated log is not created
        assert not gu.single_contig_log.exists()

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

    def test_single_contig_genome_with_multiple_contigs_raises(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "multi_contig_bin.fa.gz", ["contigX", "contigY"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path)])
        single_contig_info_path = _write_single_contig_info_tsv(
            tmp_path / "single_contig_info.tsv",
            [["ERR6769700_bin.1", "1", "Chromosome", "Circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, single_contig_info_path, test_suffix="single-contig-invalid-unittest")
        gu = GenomeUpload(args)

        with pytest.raises(ValueError, match="expected exactly 1"):
            gu.extract_genomes_info()

    def test_no_single_contig_info_provided_defaults_to_false(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("MAG1", genome_path)])
        args = _base_args(tmp_path, genome_info_path, test_suffix="no-single-contig-info-unittest")
        gu = GenomeUpload(args)
        assert gu.single_contig_metadata is None
        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        assert genome_info[alias]["single_contig"] is False
        assert genome_info[alias]["contig_id"] is None
        # the chromosome columns are always present on the genome dict, even when
        # --single-contig-info wasn't given, so compute_manifests / generate_genome_manifest
        # don't KeyError later on
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_NAME] is None
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_TYPE] is None
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_TOPOLOGY] is None
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_LOCATION] is None
        assert not gu.single_contig_log.exists()

    def test_genome_not_listed_in_single_contig_info_defaults_to_false(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("MAG1", genome_path)])
        single_contig_info_path = _write_single_contig_info_tsv(
            tmp_path / "single_contig_info.tsv",
            [["genome_that_does_not_exist", "1", "Chromosome", "Circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, single_contig_info_path, test_suffix="single-contig-unmatched-unittest")
        gu = GenomeUpload(args)
        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        # the only genome isn't listed in the single-contig info file, so it's submitted normally
        assert genome_info[alias]["single_contig"] is False
        assert genome_info[alias]["contig_id"] is None

        # the single-contig info file's genome_name doesn't match anything -> logged, not fatal
        assert gu.single_contig_log.exists()
        log_text = gu.single_contig_log.read_text()
        assert "genome_that_does_not_exist" in log_text
        assert "do not match any genome" in log_text

    def test_single_contig_genome_with_wrong_type_is_skipped_and_logged(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "single_contig_bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(
            tmp_path / "genome_info.tsv",
            [
                _genome_row("ERR6769700_bin.1", genome_path),
                _genome_row("ERR6769700_bin.2", genome_path),
            ],
        )
        single_contig_info_path = _write_single_contig_info_tsv(
            tmp_path / "single_contig_info.tsv",
            [
                ["ERR6769700_bin.1", "1", "not_a_real_type", "Circular"],
                ["ERR6769700_bin.2", "2", "Segment", "Circular"],
            ],
        )
        args = _base_args(tmp_path, genome_info_path, single_contig_info_path, test_suffix="single-contig-wrong-type-unittest")
        gu = GenomeUpload(args)

        # both genomes carry a single_contig_type that is not a SINGLE_CONTIG_CHROMOSOME_TYPES
        # key, so both are skipped and nothing is left to submit
        with pytest.raises(ValueError, match="No genomes left"):
            gu.extract_genomes_info()

        assert gu.single_contig_log.exists()
        log_text = gu.single_contig_log.read_text()
        assert "single_contig_type 'not_a_real_type' is not one of" in log_text
        assert "single_contig_type 'Segment' is not one of" in log_text
        assert "ERR6769700_bin.1" in log_text
        assert "ERR6769700_bin.2" in log_text
        assert "2 genome(s) excluded from single-contig submission" in log_text

    def test_single_contig_genome_with_wrong_name_format_is_skipped_and_logged(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "single_contig_bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path)])
        single_contig_info_path = _write_single_contig_info_tsv(
            tmp_path / "single_contig_info.tsv",
            [["ERR6769700_bin.1", "chr1", "Chromosome", "Circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, single_contig_info_path, test_suffix="single-contig-wrong-name-unittest")
        gu = GenomeUpload(args)

        # single_contig_name "chr1" is neither a digit string nor "MIT", so the only
        # genome is skipped and nothing is left to submit
        with pytest.raises(ValueError, match="No genomes left"):
            gu.extract_genomes_info()

        assert gu.single_contig_log.exists()
        log_text = gu.single_contig_log.read_text()
        assert "single_contig_name 'chr1' must be a digit string" in log_text

    def test_single_contig_mitochondrion_name_writes_chromosome_list(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "single_contig_bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path)])
        single_contig_info_path = _write_single_contig_info_tsv(
            tmp_path / "single_contig_info.tsv",
            [["ERR6769700_bin.1", "MIT", "Chromosome", "Linear"]],
        )
        args = _base_args(tmp_path, genome_info_path, single_contig_info_path, test_suffix="single-contig-mitochondrion-unittest")
        gu = GenomeUpload(args)
        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        assert genome_info[alias]["single_contig"] is True
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_NAME] == "MIT"
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_TYPE] == "Chromosome"
        assert genome_info[alias][SINGLE_CONTIG_CHROMOSOME_TOPOLOGY] == "Linear"

        genome_info[alias]["run_ref"] = "ERR6769700"
        genome_info[alias]["study"] = "ERP000000"
        genome_info[alias]["sequencingMethod"] = "ILLUMINA"

        manifest_info = compute_manifests(genome_info)[alias]
        gu.generate_genome_manifest(manifest_info, {alias: "ERS0000001"})

        chromosome_list_path = gu.manifest_dir / f"{alias}_chromosome_list.txt.gz"
        with gzip.open(chromosome_list_path, "rt") as f:
            assert f.read() == "contig1\tMIT\tLinear-Chromosome\n"

    def test_single_contig_info_missing_mandatory_column_raises(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("MAG1", genome_path)])
        # single_contig_topology is missing, even though name and type are present
        single_contig_info_path = tmp_path / "single_contig_info.tsv"
        single_contig_info_path.write_text("genome_name\tsingle_contig_name\tsingle_contig_type\nMAG1\t1\tChromosome\n")

        args = _base_args(tmp_path, genome_info_path, single_contig_info_path, test_suffix="single-contig-missing-column-unittest")
        gu = GenomeUpload(args)

        with pytest.raises(ValueError, match="single_contig_topology"):
            gu.extract_genomes_info()
