import pytest
from genomeuploader.genome_upload import GenomeUpload

@pytest.fixture
def genome_upload_instance():
    # Minimal args for instantiation; adjust as needed for your class
    args = {
        "bins": False,
        "live": False,
        "private": False,
        "tpa": False,
        "centre_name": "",
        "force": False,
        "out": ".",
        "upload_study": "",
        "genome_info": "",
        "test_suffix": None,
        "verbose": False,
    }
    return GenomeUpload(args)

@pytest.mark.parametrize("input_date,expected", [
    ("01-13-2021", "2021-01-13"),
    ("04/Jan/2010", "2010-01-04"),
])
def test_get_collection_date_reformat(genome_upload_instance, input_date, expected):
    assert genome_upload_instance.get_collection_date(input_date) == expected