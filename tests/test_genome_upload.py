import subprocess
from pathlib import Path


class Tests:

    def test_genomeuploader_end_to_end(tmp_path):
        command = [
            "python", "genomeuploader/genome_upload.py",
            "-u", "ERP159782",
            "--genome_info", "tests/fixtures/input_fixture.tsv",
            "--out", "tests/fixtures/",
            "--bins",
            "--test-unique-submission",
            "--centre_name", "EMG"
        ]

        # Run the command once
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, f"First run failed: {result.stderr}"

        # Check required output files
        expected_files = [
            "tests/fixtures/bin_upload/manifests_test/",
            "tests/fixtures/bin_upload/genome_samples.xml",
            "tests/fixtures/bin_upload/registered_bins_test.tsv",
            "tests/fixtures/bin_upload/submission.xml"
        ]
        for path in expected_files:
            assert Path(path).exists(), f"Missing expected output: {path}"

        # should have 1 line
        filepath = "tests/fixtures/bin_upload/registered_bins_test.tsv"
        with open(filepath, "r") as f:
            lines = f.readlines()
        assert len(lines) == 1