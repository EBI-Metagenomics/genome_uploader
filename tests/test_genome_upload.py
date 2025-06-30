import subprocess
from pathlib import Path
from datetime import date, datetime as dt


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
        # should have the same number of genomes
        assert len(lines) == number_of_bins
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

        with open("tests/fixtures/input_with_registrered_fixture.tsv", "r") as f:
            lines = f.readlines()
        number_of_bins2 = len(lines) - 1
        command = [
            "python",
            "genomeuploader/genome_upload.py",
            "-u",
            "ERP159782",
            "--genome_info",
            "tests/fixtures/input_with_registrered_fixture.tsv",
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
        # should have 3 line
        assert len(lines) == number_of_bins2
        # Check sample count, XML should include only new genomes, excluding registered in first round
        with open("tests/fixtures/bin_upload/genome_samples.xml") as f:
            assert f.read().count("alias=") == number_of_bins2 - number_of_bins1
