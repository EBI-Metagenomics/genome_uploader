import bz2
import gzip
from pathlib import Path

import pytest

from genomeuploader.genome_upload import *  # noqa: F401,F403 - GenomeUpload, compute_manifests
from genomeuploader.chromosome_file import ChromosomeRecord, get_fasta_contig_ids
from genomeuploader.constants import MAG_MANDATORY_FIELDS

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
    "has_chromosome",
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


def _genome_row(genome_name: str, genome_path: Path, accessions: str = "ERR6769700", has_chromosome: str = "False") -> list:
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
        has_chromosome,
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
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path, has_chromosome="True")])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["ERR6769700_bin.1", "contig1", "1", "chromosome", "circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path)
        gu = GenomeUpload(args)
        # the chromosome-info file is loaded and validated in __init__ already
        assert gu.chromosome_metadata["ERR6769700_bin.1"][0].chromosome_name == "1"

        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        assert genome_info[alias]["has_chromosome"] is True
        records = genome_info[alias]["chromosome_records"]
        assert len(records) == 1
        assert records[0].contig_name == "contig1"
        assert records[0].chromosome_name == "1"
        assert records[0].chromosome_type == "chromosome"
        assert records[0].chromosome_topology == "circular"
        assert records[0].chromosome_location == ""

        # don't call API again, force assign info for the next step
        genome_info[alias]["run_ref"] = "ERR6769700"
        genome_info[alias]["study"] = "ERP000000"
        genome_info[alias]["sequencingMethod"] = "ILLUMINA"

        manifest_info = compute_manifests(genome_info)[alias]
        gu.generate_genome_manifest(manifest_info, {alias: "ERS0000001"})

        manifest_text = (gu.manifest_dir / f"{alias}.manifest").read_text()
        chromosome_list_path = gu.manifest_dir / f"{alias}_chromosome_list.txt"

        assert chromosome_list_path.exists()
        assert f"CHROMOSOME_LIST\t{chromosome_list_path.resolve()}" in manifest_text
        assert chromosome_list_path.read_text() == "contig1\t1\tcircular-chromosome\n"

    def test_genome_with_multiple_chromosome_records_writes_all_lines(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "chromosome_bin.fa.gz", ["chr1", "chr2"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path, has_chromosome="True")])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [
                ["ERR6769700_bin.1", "chr1", "1", "chromosome", "circular"],
                ["ERR6769700_bin.1", "chr2", "2", "plasmid", "linear"],
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

        chromosome_list_path = gu.manifest_dir / f"{alias}_chromosome_list.txt"
        assert chromosome_list_path.read_text() == "chr1\t1\tcircular-chromosome\nchr2\t2\tlinear-plasmid\n"

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

    def test_genome_not_listed_in_chromosome_info_defaults_to_false(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(
            tmp_path / "genome_info.tsv",
            [_genome_row("MAG1", genome_path, has_chromosome="True"), _genome_row("MAG2", genome_path)],
        )
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["MAG1", "contig1", "1", "chromosome", "circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path, test_suffix="chromosome-unmatched-unittest")
        gu = GenomeUpload(args)
        genome_info = gu.extract_genomes_info()

        # MAG2 isn't listed in the chromosome info file, so it's submitted normally
        mag2_alias = next(alias for alias in genome_info if genome_info[alias]["genome_name"] == "MAG2")
        assert genome_info[mag2_alias]["has_chromosome"] is False
        assert genome_info[mag2_alias]["chromosome_records"] == []

    def test_chromosome_info_genome_name_not_in_genome_info_raises(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("MAG1", genome_path)])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["genome_that_does_not_exist", "contig1", "1", "chromosome", "circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path, test_suffix="chromosome-unmatched-unittest")
        gu = GenomeUpload(args)

        # a chromosome-info genome_name that doesn't match any --genome_info genome now
        # fails the whole run, rather than being silently ignored
        with pytest.raises(ValueError, match="do not match any genome") as exc_info:
            gu.extract_genomes_info()
        assert "genome_that_does_not_exist" in str(exc_info.value)

    def test_has_chromosome_true_without_chromosome_info_raises(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("MAG1", genome_path, has_chromosome="True")])
        args = _base_args(tmp_path, genome_info_path, test_suffix="has-chromosome-no-file-unittest")
        gu = GenomeUpload(args)

        with pytest.raises(ValueError, match="has_chromosome is True in --genome_info but no chromosome records") as exc_info:
            gu.extract_genomes_info()
        assert "MAG1" in str(exc_info.value)

    def test_has_chromosome_true_but_genome_not_in_chromosome_info_raises(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(
            tmp_path / "genome_info.tsv",
            [_genome_row("MAG1", genome_path, has_chromosome="True"), _genome_row("MAG2", genome_path, has_chromosome="yes")],
        )
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["MAG1", "contig1", "1", "chromosome", "circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path, test_suffix="has-chromosome-missing-records-unittest")
        gu = GenomeUpload(args)

        with pytest.raises(ValueError, match="no chromosome records were found") as exc_info:
            gu.extract_genomes_info()
        assert "for: MAG2" in str(exc_info.value)

    @pytest.mark.parametrize("has_chromosome", ["False", ""])
    def test_chromosome_records_without_has_chromosome_raises(self, tmp_path, has_chromosome):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("MAG1", genome_path, has_chromosome=has_chromosome)])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["MAG1", "contig1", "1", "chromosome", "circular"]],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path, test_suffix="records-without-flag-unittest")
        gu = GenomeUpload(args)

        with pytest.raises(ValueError, match="has_chromosome is not True in --genome_info for: MAG1"):
            gu.extract_genomes_info()

    def test_has_chromosome_column_missing_defaults_to_false(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        row = _genome_row("MAG1", genome_path)[:-1]
        genome_info_path = _write_tsv(tmp_path / "genome_info.tsv", GENOME_TSV_HEADER[:-1], [row])
        args = _base_args(tmp_path, genome_info_path, test_suffix="has-chromosome-missing-column-unittest")
        genome_info = GenomeUpload(args).extract_genomes_info()

        alias = next(iter(genome_info))
        assert genome_info[alias]["has_chromosome"] is False

    def test_mags_empty_has_chromosome_is_not_reported_as_missing(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("MAG1", genome_path, has_chromosome="")])
        args = _base_args(tmp_path, genome_info_path, test_suffix="mags-empty-has-chromosome-unittest")
        args["bins"] = False
        mandatory_fields_before = list(MAG_MANDATORY_FIELDS)

        genome_info = GenomeUpload(args).validate_metadata_tsv()

        # has_chromosome is optional, and reading it must not add it to the shared constant
        assert next(iter(genome_info.values()))["has_chromosome"] is False
        assert MAG_MANDATORY_FIELDS == mandatory_fields_before

    def test_has_chromosome_invalid_value_raises(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("MAG1", genome_path, has_chromosome="maybe")])
        args = _base_args(tmp_path, genome_info_path, test_suffix="has-chromosome-invalid-unittest")

        with pytest.raises(ValueError, match="has_chromosome must be one of"):
            GenomeUpload(args).extract_genomes_info()

    def test_bool_columns_accept_true_and_false_values_case_insensitively(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        row = _genome_row("MAG1", genome_path)
        row[GENOME_TSV_HEADER.index("co-assembly")] = "no"
        row[GENOME_TSV_HEADER.index("rRNA_presence")] = " YES "
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [row])
        args = _base_args(tmp_path, genome_info_path, test_suffix="bool-values-unittest")
        genome_info = GenomeUpload(args).validate_metadata_tsv()

        values = next(iter(genome_info.values()))
        assert values["co-assembly"] is False
        assert values["rRNA_presence"] is True
        assert values["has_chromosome"] is False

    @pytest.mark.parametrize("column", ["co-assembly", "rRNA_presence"])
    def test_bool_column_invalid_value_raises(self, tmp_path, column):
        genome_path = _write_fasta_gz(tmp_path / "bin.fa.gz", ["contig1"])
        row = _genome_row("MAG1", genome_path)
        row[GENOME_TSV_HEADER.index(column)] = "maybe"
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [row])
        args = _base_args(tmp_path, genome_info_path, test_suffix="bool-invalid-unittest")

        with pytest.raises(ValueError, match=f"{column} must be one of .* for the following genomes: MAG1"):
            GenomeUpload(args).validate_metadata_tsv()

    def test_chromosome_info_contig_name_not_in_fasta_raises(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "chromosome_bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path, has_chromosome="True")])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["ERR6769700_bin.1", "contig_that_does_not_exist", "1", "chromosome", "circular"]],
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
                ["ERR6769700_bin.1", "contig1", "1", "chromosome", "circular"],
                ["ERR6769700_bin.1", "contig1", "2", "plasmid", "linear"],
            ],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-duplicate-contig-unittest")

        with pytest.raises(ValueError, match="is listed in more than one chromosome record"):
            GenomeUpload(args)

    def test_chromosome_info_with_wrong_type_raises(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [
                ["ERR6769700_bin.1", "contig1", "1", "not_a_real_type", "circular"],
                ["ERR6769700_bin.2", "contig1", "2", "Segment", "circular"],
            ],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-wrong-type-unittest")

        # both rows carry a chromosome_type that is not one of CHROMOSOME_TYPES, so the
        # whole file is rejected before anything else (--genome_info doesn't even need to exist)
        with pytest.raises(ValueError, match="failed validation") as exc_info:
            GenomeUpload(args)

        assert "chromosome_type 'not_a_real_type' is not one of" in str(exc_info.value)
        assert "chromosome_type 'segment' is not one of" in str(exc_info.value)
        assert "ERR6769700_bin.1" in str(exc_info.value)
        assert "ERR6769700_bin.2" in str(exc_info.value)

    def test_chromosome_info_with_wrong_name_format_raises(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["ERR6769700_bin.1", "contig1", "chr1", "chromosome", "circular"]],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-wrong-name-unittest")

        # chromosome_name "chr1" contains the forbidden substring "chr"
        with pytest.raises(ValueError, match="must not contain 'chr'"):
            GenomeUpload(args)

    def test_chromosome_info_chromosome_name_rules_raise(self, tmp_path):
        long_name = "a" * 33
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [
                ["ERR6769700_bin.1", "contig1", "_1", "chromosome", "circular"],
                ["ERR6769700_bin.1", "contig2", "1$", "chromosome", "circular"],
                ["ERR6769700_bin.1", "contig3", long_name, "chromosome", "circular"],
                ["ERR6769700_bin.1", "contig4", "pPLASMID1", "plasmid", "circular"],
                ["ERR6769700_bin.1", "contig5", "Linkage-Group2", "linkage_group", "linear"],
            ],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-name-rules-unittest")

        with pytest.raises(ValueError, match="failed validation") as exc_info:
            GenomeUpload(args)

        message = str(exc_info.value)
        assert "chromosome_name '_1' must match the pattern" in message
        assert "chromosome_name '1$' must match the pattern" in message
        assert f"chromosome_name '{long_name}' must be shorter than 33 characters" in message
        assert "chromosome_name 'pPLASMID1' must not contain 'plasmid'" in message
        assert "chromosome_name 'Linkage-Group2' must not contain 'linkage-group'" in message

    def test_chromosome_info_valid_chromosome_names_load(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [
                ["ERR6769700_bin.1", "contig1", "1", "chromosome", "circular"],
                ["ERR6769700_bin.1", "contig2", "MIT", "chromosome", "linear"],
                ["ERR6769700_bin.1", "contig3", "pA_2#b-c.d", "plasmid", "circular"],
                ["ERR6769700_bin.1", "contig4", "a" * 32, "plasmid", "circular"],
            ],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-valid-names-unittest")

        gu = GenomeUpload(args)
        names = [r.chromosome_name for r in gu.chromosome_metadata["ERR6769700_bin.1"]]
        assert names == ["1", "MIT", "pA_2#b-c.d", "a" * 32]

    def test_chromosome_info_duplicate_chromosome_name_for_genome_raises(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [
                ["ERR6769700_bin.1", "contig1", "1", "chromosome", "circular"],
                ["ERR6769700_bin.1", "contig2", "1", "plasmid", "circular"],
                # the same name in a different genome (assembly) is fine
                ["ERR6769700_bin.2", "contig1", "1", "chromosome", "circular"],
            ],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-duplicate-name-unittest")

        with pytest.raises(ValueError, match="failed validation") as exc_info:
            GenomeUpload(args)

        message = str(exc_info.value)
        assert "Genome 'ERR6769700_bin.1': chromosome_name '1' is listed in more than one chromosome record" in message
        assert "ERR6769700_bin.2" not in message

    def test_chromosome_type_and_topology_are_case_insensitive(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [
                ["ERR6769700_bin.1", "contig1", "1", "CHROMOSOME", "Circular"],
                ["ERR6769700_bin.1", "contig2", "2", "Linkage_Group", "LINEAR"],
            ],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-case-insensitive-unittest")

        gu = GenomeUpload(args)
        records = gu.chromosome_metadata["ERR6769700_bin.1"]
        assert [r.chromosome_type for r in records] == ["chromosome", "linkage_group"]
        assert [r.chromosome_topology for r in records] == ["circular", "linear"]

    def test_chromosome_info_empty_contig_name_raises(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["ERR6769700_bin.1", "", "1", "chromosome", "circular"]],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-empty-contig-unittest")

        with pytest.raises(ValueError, match="contig_name must not be empty"):
            GenomeUpload(args)

    def test_chromosome_info_empty_mandatory_fields_raises(self, tmp_path):
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [
                ["", "contig1", "1", "chromosome", "circular"],
                ["ERR6769700_bin.1", "contig2", "", "", "linear"],
            ],
        )
        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-empty-fields-unittest")

        with pytest.raises(ValueError, match="failed validation") as exc_info:
            GenomeUpload(args)

        message = str(exc_info.value)
        assert "Row 2 (genome ''): mandatory field(s) genome_name must not be empty" in message
        assert "Row 3 (genome 'ERR6769700_bin.1'): mandatory field(s) chromosome_name, chromosome_type must not be empty" in message
        # empty values are reported once, not again as vocabulary mismatches
        assert "is not one of" not in message
        assert "must match the pattern" not in message

    def test_mitochondrion_name_writes_chromosome_list(self, tmp_path):
        genome_path = _write_fasta_gz(tmp_path / "chromosome_bin.fa.gz", ["contig1"])
        genome_info_path = _write_genome_tsv(tmp_path / "genome_info.tsv", [_genome_row("ERR6769700_bin.1", genome_path, has_chromosome="True")])
        chromosome_info_path = _write_chromosome_info_tsv(
            tmp_path / "chromosome_info.tsv",
            [["ERR6769700_bin.1", "contig1", "MIT", "chromosome", "linear"]],
        )
        args = _base_args(tmp_path, genome_info_path, chromosome_info_path, test_suffix="chromosome-mitochondrion-unittest")
        gu = GenomeUpload(args)
        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        records = genome_info[alias]["chromosome_records"]
        assert genome_info[alias]["has_chromosome"] is True
        assert records[0].chromosome_name == "MIT"
        assert records[0].chromosome_type == "chromosome"
        assert records[0].chromosome_topology == "linear"

        genome_info[alias]["run_ref"] = "ERR6769700"
        genome_info[alias]["study"] = "ERP000000"
        genome_info[alias]["sequencingMethod"] = "ILLUMINA"

        manifest_info = compute_manifests(genome_info)[alias]
        gu.generate_genome_manifest(manifest_info, {alias: "ERS0000001"})

        chromosome_list_path = gu.manifest_dir / f"{alias}_chromosome_list.txt"
        assert chromosome_list_path.read_text() == "contig1\tMIT\tlinear-chromosome\n"

    def test_chromosome_info_missing_mandatory_column_raises(self, tmp_path):
        # chromosome_topology is missing, even though the others are present
        chromosome_info_path = tmp_path / "chromosome_info.tsv"
        chromosome_info_path.write_text(
            "genome_name\tcontig_name\tchromosome_name\tchromosome_type\nMAG1\tcontig1\t1\tchromosome\n"
        )

        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-missing-column-unittest")

        with pytest.raises(ValueError, match="chromosome_topology"):
            GenomeUpload(args)

    def test_chromosome_info_missing_contig_name_column_raises(self, tmp_path):
        chromosome_info_path = tmp_path / "chromosome_info.tsv"
        chromosome_info_path.write_text(
            "genome_name\tchromosome_name\tchromosome_type\tchromosome_topology\nMAG1\t1\tchromosome\tcircular\n"
        )

        args = _base_args(tmp_path, tmp_path / "genome_info.tsv", chromosome_info_path, test_suffix="chromosome-missing-contig-column-unittest")

        with pytest.raises(ValueError, match="contig_name"):
            GenomeUpload(args)


@pytest.mark.parametrize("suffix, opener", [(".gz", gzip.open), (".bz2", bz2.open)])
def test_get_fasta_contig_ids_reads_compressed_fasta(tmp_path, suffix, opener):
    fasta_path = tmp_path / f"genome.fa{suffix}"
    with opener(fasta_path, "wt") as f:
        f.write(">contig1\nACGT\n>contig2\nGGCC\n")

    assert get_fasta_contig_ids(fasta_path) == ["contig1", "contig2"]


def test_get_fasta_contig_ids_rejects_uncompressed_extension(tmp_path):
    fasta_path = tmp_path / "genome.fa"
    fasta_path.write_text(">contig1\nACGT\n")

    with pytest.raises(ValueError, match="must be compressed with gzip or bzip2"):
        get_fasta_contig_ids(fasta_path)


@pytest.mark.parametrize("suffix", [".gz", ".bz2"])
def test_get_fasta_contig_ids_rejects_invalid_compressed_data(tmp_path, suffix):
    fasta_path = tmp_path / f"genome.fa{suffix}"
    fasta_path.write_text(">contig1\nACGT\n")

    with pytest.raises(ValueError, match="is not valid"):
        get_fasta_contig_ids(fasta_path)


def _record(**overrides) -> ChromosomeRecord:
    row = {
        "genome_name": "MAG1",
        "contig_name": "contig1",
        "chromosome_name": "1",
        "chromosome_type": "chromosome",
        "chromosome_topology": "circular",
        "chromosome_location": "",
    }
    row.update(overrides)
    return ChromosomeRecord.from_row(row)


def test_chromosome_record_from_row_lowercases_type_and_topology():
    record = _record(chromosome_type="Linkage_Group", chromosome_topology="LINEAR")
    assert record.chromosome_type == "linkage_group"
    assert record.chromosome_topology == "linear"
    assert record.validate() == []


def test_chromosome_record_is_immutable():
    record = _record()
    with pytest.raises(AttributeError):
        record.chromosome_name = "2"


@pytest.mark.parametrize(
    "name, expected",
    [
        ("1", None),
        ("MIT", None),
        ("pA_2#b-c.d", None),
        ("", None),
        ("_1", "must match the pattern"),
        ("1$", "must match the pattern"),
        ("a" * 33, "must be shorter than 33 characters"),
        ("Chr1", "must not contain 'chr'"),
        ("pPLASMID", "must not contain 'plasmid'"),
    ],
)
def test_chromosome_record_validate_chromosome_name(name, expected):
    errors = _record(chromosome_name=name).validate_chromosome_name()
    if expected is None:
        assert errors == []
    else:
        assert any(expected in error for error in errors)


@pytest.mark.parametrize("value, valid", [("plasmid", True), ("Segmented", True), ("", True), ("segment", False)])
def test_chromosome_record_validate_chromosome_type(value, valid):
    assert (_record(chromosome_type=value).validate_chromosome_type() == []) is valid


@pytest.mark.parametrize("value, valid", [("linear", True), ("Circular", True), ("", True), ("loop", False)])
def test_chromosome_record_validate_chromosome_topology(value, valid):
    assert (_record(chromosome_topology=value).validate_chromosome_topology() == []) is valid


@pytest.mark.parametrize("value, valid", [("Mitochondrion", True), ("", True), ("mitochondrion", False)])
def test_chromosome_record_validate_chromosome_location(value, valid):
    assert (_record(chromosome_location=value).validate_chromosome_location() == []) is valid


def test_chromosome_record_to_chromosome_list_line():
    assert _record().to_chromosome_list_line() == "contig1\t1\tcircular-chromosome"
    assert (
        _record(chromosome_name="MIT", chromosome_topology="linear", chromosome_location="Mitochondrion").to_chromosome_list_line()
        == "contig1\tMIT\tlinear-chromosome\tMitochondrion"
    )
