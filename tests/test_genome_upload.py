import subprocess
from datetime import datetime as dt
from pathlib import Path

import pytest
import responses as responses_lib

from genomeuploader.ena_submit import EnaSubmit
from genomeuploader.genome_upload import *


class Tests:
    @responses_lib.activate
    def test_register_genome_samples_in_ena_uses_webin_v2_queue(self, tmp_path, monkeypatch):
        monkeypatch.setattr("genomeuploader.ena_submit.CredentialsManager.get_credentials", lambda: ("user", "password"))
        monkeypatch.setattr("genomeuploader.ena_submit.time.sleep", lambda _: None)

        payload_xml = tmp_path / "submission.xml"
        payload_xml.write_text(
            "<WEBIN><SUBMISSION_SET><SUBMISSION><ACTIONS><ACTION><ADD/></ACTION></ACTIONS></SUBMISSION></SUBMISSION_SET><SAMPLE_SET><SAMPLE alias=\"genome_1\" center_name=\"EMG\"/></SAMPLE_SET></WEBIN>"
        )

        receipt_path = tmp_path / "submission_receipt.xml"
        queue_url = "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit/queue"
        poll_url = "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit/poll/ERA000001"
        final_receipt = '<RECEIPT success="true"><SAMPLE accession="ERS000001" alias="genome_1"/></RECEIPT>'

        responses_lib.add(responses_lib.POST, queue_url, json={"submissionId": "ERA000001", "_links": {"poll": {"href": poll_url}}}, status=200)
        responses_lib.add(responses_lib.GET, poll_url, body="queued", status=202)
        responses_lib.add(responses_lib.GET, poll_url, body=final_receipt, status=200, content_type="application/xml")

        ena_submit = EnaSubmit(payload_xml, receipt_path, 1, live=False)
        alias_map = ena_submit.register_genome_samples_in_ena()

        assert alias_map == {"genome_1": "ERS000001"}
        assert receipt_path.read_text() == final_receipt
        assert b"<WEBIN>" in responses_lib.calls[0].request.body


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
            "tests/fixtures/bin_upload/submission_receipt.xml",
            "tests/fixtures/bin_upload/registered_bins_test.tsv",
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

    def test_genomeuploader_re_register_samples(tmp_path):
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
            "tests/fixtures/bin_upload/genome_samples_retry.xml",
            "tests/fixtures/bin_upload/submission_receipt.xml",
            "tests/fixtures/bin_upload/submission_receipt_retry.xml",
            "tests/fixtures/bin_upload/registered_bins_test.tsv",
        ]
        for path in expected_files:
            assert Path(path).exists(), f"Missing expected output: {path}"

        # check registered samples tsv
        filepath = "tests/fixtures/bin_upload/registered_bins_test.tsv"
        with open(filepath, "r") as f:
            lines = f.readlines()
        # should have 3 lines + header
        assert len(lines) == number_of_bins2 + 1
        # Check sample count, retry XML should include only new genomes, excluding registered in first round
        # genome_samples will still contain all of them
        with open("tests/fixtures/bin_upload/genome_samples.xml") as f:
            assert f.read().count("alias=") == number_of_bins2
        with open("tests/fixtures/bin_upload/genome_samples_retry.xml") as f:
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
        responses_lib.add(
            responses_lib.GET,
            "https://www.ebi.ac.uk/ena/browser/api/xml/ERZ23498047",
            body='<ANALYSIS_SET><ANALYSIS><RUN_REF accession="ERR4918394"/></ANALYSIS></ANALYSIS_SET>',
            content_type="application/xml",
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
        assert genome_info["MAG1"]["run_ref"] == "ERR4918394"

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
        assert genome_info["MAG1"]["run_ref"] == "SRR12059190"

    def test_generate_genome_manifest_excludes_run_ref_when_missing(self, tmp_path):
        args = {
            "bins": True,
            "live": False,
            "private": False,
            "tpa": False,
            "centre_name": "EMG",
            "force": False,
            "out": str(tmp_path),
            "upload_study": "ERP000000",
            "genome_info": str(tmp_path / "genome_info.tsv"),
            "test_suffix": None,
        }
        gu = GenomeUpload(args)
        genome_info = {
            "genome_name": "MAG1",
            "alias": "MAG1",
            "coverageDepth": 10.0,
            "assembler": "megahit_v1.0.0",
            "sequencingMethod": "ILLUMINA",
            "study": "PRJEB000000",
            "accessionType": "assembly",
            "accessions": "ERZ0000001",
            "genome_path": str(tmp_path / "MAG1.fa.gz"),
            "co-assembly": False,
            "run_ref": None,
        }
        alias_to_sample = {"MAG1": "ERS0000001"}

        gu.generate_genome_manifest(genome_info, alias_to_sample)
        manifest_text = (gu.manifest_dir / "MAG1.manifest").read_text()

        assert "RUN_REF\t" not in manifest_text

    def test_single_contig_genome_extracts_contig_id_and_writes_chromosome_list(self, tmp_path):
        args = {
            "bins": True,
            "live": False,
            "private": False,
            "tpa": False,
            "centre_name": "EMG",
            "force": False,
            "out": str(tmp_path),
            "upload_study": "ERP000000",
            "genome_info": "tests/fixtures/input_single_contig_fixture.tsv",
            "test_suffix": "single-contig-unittest",
        }
        gu = GenomeUpload(args)
        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        assert genome_info[alias]["single_contig"] is True
        assert genome_info[alias]["contig_id"] == "contig1"

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
        args = {
            "bins": True,
            "live": False,
            "private": False,
            "tpa": False,
            "centre_name": "EMG",
            "force": False,
            "out": str(tmp_path),
            "upload_study": "ERP000000",
            "genome_info": "tests/fixtures/input_single_contig_invalid_fixture.tsv",
            "test_suffix": "single-contig-invalid-unittest",
        }
        gu = GenomeUpload(args)

        with pytest.raises(ValueError, match="expected exactly 1"):
            gu.extract_genomes_info()

    def test_single_contig_column_absent_defaults_to_false(self, tmp_path):
        args = {
            "bins": True,
            "live": False,
            "private": False,
            "tpa": False,
            "centre_name": "EMG",
            "force": False,
            "out": str(tmp_path),
            "upload_study": "ERP000000",
            "genome_info": "tests/fixtures/input_fixture.tsv",
            "test_suffix": "no-single-contig-column-unittest",
        }
        gu = GenomeUpload(args)
        genome_info = gu.extract_genomes_info()
        alias = next(iter(genome_info))

        assert genome_info[alias]["single_contig"] is False
        assert genome_info[alias]["contig_id"] is None