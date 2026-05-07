import json

import pytest
import responses

from genomeuploader.ena import EnaQuery
from genomeuploader.genome_upload import GenomeUpload
from test_ena import read_json

@pytest.fixture
def genome_upload_public_instance(tmp_path):
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
        "verbose": False,
    }
    return GenomeUpload(args)

@pytest.fixture
def genome_upload_private_instance(tmp_path):
    args = {
        "bins": True,
        "live": False,
        "private": True,
        "tpa": False,
        "centre_name": "EMG",
        "force": False,
        "out": str(tmp_path),
        "upload_study": "ERP000000",
        "genome_info": str(tmp_path / "genome_info.tsv"),
        "test_suffix": None,
        "verbose": False,
    }
    return GenomeUpload(args)

@responses.activate
def test_all_public_get_metadata_with_privacy_success(genome_upload_public_instance):
    sample_accession = "SAMN00000001"
    responses.add(
        responses.POST,
        "https://www.ebi.ac.uk/ena/portal/api/search",
        json=[{
            "sample_accession": sample_accession,
            "collection_date": "2022-01-01",
            "country": "France",
            "location": "48.8566 N 2.3522 E"
        }]
    )

    metadata = genome_upload_public_instance.get_metadata_with_privacy(sample_accession, privacy_status=False)
    assert metadata["country"] == "France"
    assert metadata["latitude"].startswith("48.85")
    assert metadata["longitude"].startswith("2.35")
    assert metadata["collection_date"] == "2022-01-01"

@responses.activate
def test_get_metadata_with_privacy_missing_fields(genome_upload_public_instance):
    # Mock ENA sample API response with missing fields
    sample_accession = "SAMN00000002"
    responses.add(
        responses.POST,
        "https://www.ebi.ac.uk/ena/portal/api/search",
        json=[{
            "sample_accession": sample_accession,
            "collection_date": "",
            "country": "",
            "location": ""
        }]
    )

    metadata = genome_upload_public_instance.get_metadata_with_privacy(sample_accession, privacy_status=False)

    assert metadata["country"] == "missing: third party data"
    assert metadata["latitude"] == "missing: third party data"
    assert metadata["longitude"] == "missing: third party data"
    assert metadata["collection_date"] == "missing: third party data"

@responses.activate
def test_get_sample_metadata_private_submission_public_sample(monkeypatch, genome_upload_private_instance, public_sample_data):
    # Simulate first privacy attempt returns None, second returns valid data
    def emulate_get_metadata_with_privacy(accession, privacy_status):
        if privacy_status is False:
            return None
        else:
            metadata_dict = genome_upload_private_instance.get_location_metadata(public_sample_data)
            if metadata_dict:
                metadata_dict["collection_date"] = genome_upload_private_instance.get_collection_date(public_sample_data["collection_date"])
                return metadata_dict
            else:
                return None

    monkeypatch.setattr(genome_upload_private_instance, "get_metadata_with_privacy", emulate_get_metadata_with_privacy)
    metadata = genome_upload_private_instance.get_sample_metadata("SAMEA7687881")

    assert metadata["country"] == "Norway"
    assert metadata["latitude"].startswith("66.07")
    assert metadata["longitude"].startswith("12.58")
    assert metadata["collection_date"] == "2019-06-11"