import subprocess
from datetime import datetime as dt
from pathlib import Path

import responses as responses_lib

from genomeuploader.genome_upload import *


class Tests:
    def test_genomeuploader_end_to_end(tmp_path):
        timestamp = str(int(dt.timestamp(dt.now())))
        with open("tests/fixtures/input_fixture.tsv", "r") as f:
            lines = f.readlines()
        number_of_bins = len(lines) - 1
        command = [
            "python",
            "genomeuploader/genome_upload.py",
            "-u",
            "ERP159782",
            "--genome_info",
            "tests/fixtures/input_fixture.tsv",
            "--out",
            "tests/fixtures/",
            "--bins",
            "--test-suffix",
            f"end-to-end-{timestamp}",
            "--centre_name",
            "EMG",
        ]

        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"Run failed: {result.stderr}"

        # Check required output files
        expected_files = [
            "tests/fixtures/bin_upload/manifests_test/",
            "tests/fixtures/bin_upload/genome_samples.xml",
            "tests/fixtures/bin_upload/registered_bins_test.tsv",
            "tests/fixtures/bin_upload/submission.xml",
        ]
        for path in expected_files:
            assert Path(path).exists(), f"Missing expected output: {path}"

        # check registered samples tsv
        filepath = "tests/fixtures/bin_upload/registered_bins_test.tsv"
        with open(filepath, "r") as f:
            lines = f.readlines()
        # should have the same number of genomes (including header)
        assert len(lines) == number_of_bins + 1
        # should have sample id (ERS) and suffix from --test-suffix command
        assert "ERS" in "".join(lines) and "end-to-end" in "".join(lines)

    def test_genomeuploader_registered_samples(tmp_path):
        timestamp = str(int(dt.timestamp(dt.now())))
        with open("tests/fixtures/input_fixture.tsv", "r") as f:
            lines = f.readlines()
        number_of_bins1 = len(lines) - 1
        command = [
            "python",
            "genomeuploader/genome_upload.py",
            "-u",
            "ERP159782",
            "--genome_info",
            "tests/fixtures/input_fixture.tsv",
            "--out",
            "tests/fixtures/",
            "--bins",
            "--test-suffix",
            f"registered-{timestamp}",
            "--centre_name",
            "EMG",
        ]
        result1 = subprocess.run(command, capture_output=True, text=True)
        assert result1.returncode == 0, f"First run failed: {result1.stderr}"

        with open("tests/fixtures/input_with_registered_fixture.tsv", "r") as f:
            lines = f.readlines()
        number_of_bins2 = len(lines) - 1
        command = [
            "python",
            "genomeuploader/genome_upload.py",
            "-u",
            "ERP159782",
            "--genome_info",
            "tests/fixtures/input_with_registered_fixture.tsv",
            "--out",
            "tests/fixtures/",
            "--bins",
            "--test-suffix",
            f"registered-{timestamp}",
            "--centre_name",
            "EMG",
        ]
        result2 = subprocess.run(command, capture_output=True, text=True)
        assert result2.returncode == 0, f"Second run failed: {result2.stderr}"

        # Check required output files
        expected_files = [
            "tests/fixtures/bin_upload/manifests_test/",
            "tests/fixtures/bin_upload/genome_samples.xml",
            "tests/fixtures/bin_upload/registered_bins_test.tsv",
            "tests/fixtures/bin_upload/submission.xml",
        ]
        for path in expected_files:
            assert Path(path).exists(), f"Missing expected output: {path}"

        # check registered samples tsv
        filepath = "tests/fixtures/bin_upload/registered_bins_test.tsv"
        with open(filepath, "r") as f:
            lines = f.readlines()
        # should have 3 lines + header
        assert len(lines) == number_of_bins2 + 1
        # Check sample count, XML should include only new genomes, excluding registered in first round
        with open("tests/fixtures/bin_upload/genome_samples.xml") as f:
            assert f.read().count("alias=") == number_of_bins2 - number_of_bins1

    @responses_lib.activate
    def test_genomeuploader_extract_ena_info_assembly(tmp_path):
        genome_info = {
            "MAG1": {
                "accessions": ["ERZ23498047"],
                "accessionType": "assembly",
                "co-assembly": False,
            }
        }

        # 1. assembly_info query
        responses_lib.add(
            responses_lib.POST,
            "https://www.ebi.ac.uk/ena/portal/api/search",
            json=[{"study_accession": "PRJEB71644", "sample_accession": "SAMEA114545946", "sampling_platform": "ILLUMINA"}],
        )
        # 2. study query (called by get_project_description)
        responses_lib.add(
            responses_lib.POST,
            "https://www.ebi.ac.uk/ena/portal/api/search",
            json=[{"study_accession": "PRJEB71644", "study_description": "Test metagenome study"}],
        )
        # 3. sample query (called by get_sample_metadata)
        responses_lib.add(
            responses_lib.POST,
            "https://www.ebi.ac.uk/ena/portal/api/search",
            json=[{"sample_accession": "SAMEA114545946", "collection_date": "2021-01-01", "country": "USA", "location": "40.7128 N 74.0060 W"}],
        )

        args = {
            "bins": True,
            "live": False,
            "private": False,
            "tpa": False,
            "centre_name": "EMG",
            "force": False,
            "out": str(tmp_path),
            "upload_study": "ERP000000",
            "genome_info": str(tmp_path) + "genome_info.tsv",
            "test_suffix": None,
        }
        gu = GenomeUpload(args)
        gu.extract_ena_info(genome_info)

        assert genome_info["MAG1"]["study"] == "PRJEB71644"
        assert genome_info["MAG1"]["sequencingMethod"] == "ILLUMINA"
        assert genome_info["MAG1"]["description"] == "Test metagenome study"
        assert genome_info["MAG1"]["sample_accessions"] == "SAMEA114545946"
        assert genome_info["MAG1"]["collectionDate"] == "2021-01-01"
        assert genome_info["MAG1"]["country"] == "USA"
        assert genome_info["MAG1"]["accessions"] == "ERZ23498047"

    @responses_lib.activate
    def test_genomeuploader_extract_ena_info_run(tmp_path):
        genome_info = {
            "MAG1": {
                "accessions": ["SRR12059190"],
                "accessionType": "run",
                "co-assembly": False,
            }
        }

        # 1. run query
        responses_lib.add(
            responses_lib.POST,
            "https://www.ebi.ac.uk/ena/portal/api/search",
            json=[{"run_accession": "SRR12059190", "secondary_study_accession": "SRP272267", "sample_accession": "SAMN15965141"}],
        )
        # 2. study query (called by get_project_description)
        responses_lib.add(
            responses_lib.POST,
            "https://www.ebi.ac.uk/ena/portal/api/search",
            json=[{"study_accession": "SRP272267", "study_description": "Test run study"}],
        )
        # 3. sample query (called by get_sample_metadata)
        responses_lib.add(
            responses_lib.POST,
            "https://www.ebi.ac.uk/ena/portal/api/search",
            json=[{"sample_accession": "SAMN15965141", "collection_date": "2020-06-15", "country": "Germany", "location": "52.5200 N 13.4050 E"}],
        )

        args = {
            "bins": True,
            "live": False,
            "private": False,
            "tpa": False,
            "centre_name": "EMG",
            "force": False,
            "out": str(tmp_path),
            "upload_study": "ERP000000",
            "genome_info": str(tmp_path) + "genome_info.tsv",
            "test_suffix": None,
        }
        gu = GenomeUpload(args)
        gu.extract_ena_info(genome_info)

        assert genome_info["MAG1"]["study"] == "SRP272267"
        assert genome_info["MAG1"]["sequencingMethod"] is None
        assert genome_info["MAG1"]["description"] == "Test run study"
        assert genome_info["MAG1"]["sample_accessions"] == "SAMN15965141"
        assert genome_info["MAG1"]["collectionDate"] == "2020-06-15"
        assert genome_info["MAG1"]["country"] == "Germany"
        assert genome_info["MAG1"]["accessions"] == "SRR12059190"