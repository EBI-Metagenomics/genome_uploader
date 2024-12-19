from pathlib import Path

import pytest

@pytest.fixture(scope="module")
def test_file_dir():
    return Path(__file__).resolve().parent / Path("fixtures/responses")

@pytest.fixture(scope="module")
def public_study_json(test_file_dir):
    return test_file_dir / Path("public_study.json")

@pytest.fixture(scope="module")
def private_study_xml(test_file_dir):
    return test_file_dir / Path("private_study.xml")

@pytest.fixture(scope="module")
def public_run_json(test_file_dir):
    return test_file_dir / Path("public_run.json")

@pytest.fixture(scope="module")
def private_run_json(test_file_dir):
    return test_file_dir / Path("private_run.json")

@pytest.fixture(scope="module")
def public_run_from_assembly_xml(test_file_dir):
    return test_file_dir / Path("public_run_from_assembly.xml")

@pytest.fixture(scope="module")
def private_run_from_assembly_xml(test_file_dir):
    return test_file_dir / Path("private_run_from_assembly.xml")

@pytest.fixture(scope="module")
def public_study_runs_json(test_file_dir):
    return test_file_dir / Path("public_study_runs.json")

@pytest.fixture(scope="module")
def private_study_runs_json(test_file_dir):
    return test_file_dir / Path("private_study_runs.json")

@pytest.fixture(scope="module")
def public_sample_json(test_file_dir):
    return test_file_dir / Path("public_sample.json")

@pytest.fixture(scope="module")
def private_sample_xml(test_file_dir):
    return test_file_dir / Path("private_sample.xml")


@pytest.fixture
def public_study_data():
    return {
    "study_accession": "PRJEB41657",
    "study_description": "Metagenomic raw reads, assemblies, and bins derived from HoloFood salmon gut samples from trial A and B. The samples in this project contributed to the salmon MAG catalogue (project: PRJEB55376 [ERP140265])."
    }

#   private XML has no PRJ study ID
@pytest.fixture
def private_study_data():
    return {
    "study_accession": "ERP125469",
    "study_description": "Metagenomic raw reads, assemblies, and bins derived from HoloFood salmon gut samples from trial A and B. The samples in this project contributed to the salmon MAG catalogue (project: PRJEB55376 [ERP140265])."
    }

@pytest.fixture
def public_run_data():
    return {
        "run_accession": "ERR4918394",
        "sample_accession": "SAMEA7687881",
        "secondary_study_accession": "ERP125469"
    }

#   private json only has ERS sample IDs
@pytest.fixture
def private_run_data():
    return {
        "run_accession": "ERR4918394",
        "sample_accession": "ERS5444411",
        "secondary_study_accession": "ERP125469"
    }

@pytest.fixture
def public_sample_data():
    return {
        "sample_accession": "SAMEA7687881",
        "collection_date": "2019-06-11",
        "country": "Norway",
        "location": "66.079905 N 12.587848 E"
  }

#  private xml can only query with RS accession and lat and long are split
@pytest.fixture
def private_sample_data():
    return {
        "sample_accession": "ERS5444411",
        "collection_date": "2019-06-11",
        "country": "Norway",
        'latitude': '66.079905',
        'longitude': '12.587848'
  }
