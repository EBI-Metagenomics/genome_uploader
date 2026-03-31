from genomeuploader.genome_upload import GenomeUpload


def _mock_genome(alias: str) -> dict:
    return {
        "alias": alias,
        "co-assembly": False,
        "accessions": "ERR000001",
        "study": "ERP000001",
        "taxID": "2",
        "scientific_name": "Bacteria",
        "description": "test project",
        "sequencingMethod": "Illumina",
        "assembly_software": "megahit_v1.2.9",
        "MAG_quality": "Multiple fragments",
        "binning_software": "metabat2_v2.15",
        "binning_parameters": "default",
        "stats_generation_software": "checkm2_v1.0",
        "completeness": 90.0,
        "contamination": 2.0,
        "isolationSource": "marine metagenome",
        "collectionDate": "2024-01-01",
        "country": "Norway",
        "latitude": "66.079905",
        "longitude": "12.587848",
        "broad_environment": "marine biome",
        "local_environment": "coastal water",
        "environmental_medium": "water",
        "sample_accessions": "SAMEA000001",
        "metagenome": "marine metagenome",
    }


def test_split_genomes_by_sample_xml_size(tmp_path):
    args = {
        "upload_study": "ERP000001",
        "genome_info": "tests/fixtures/input_fixture.tsv",
        "mags": False,
        "bins": True,
        "out": str(tmp_path),
        "force": False,
        "live": False,
        "test_suffix": "unit",
        "tpa": False,
        "centre_name": "EMG",
        "private": False,
    }

    uploader = GenomeUpload(args)
    uploader.max_sample_xml_size_bytes = 4300

    genomes = {
        "genome_1": _mock_genome("genome_1"),
        "genome_2": _mock_genome("genome_2"),
        "genome_3": _mock_genome("genome_3"),
        "genome_4": _mock_genome("genome_4"),
    }

    batches = uploader.split_genomes_by_sample_xml_size(genomes)

    assert len(batches) > 1
    assert sum(len(batch_aliases) for _, batch_aliases in batches) == len(genomes)
    assert not uploader.samples_xml.exists()

    for batch_xml_bytes, batch_info in batches:
        assert isinstance(batch_xml_bytes, bytes)
        assert len(batch_info) > 0
        assert len(batch_xml_bytes) < uploader.max_sample_xml_size_bytes


def test_write_genomes_xml_paths(tmp_path):
    args = {
        "upload_study": "ERP000001",
        "genome_info": "tests/fixtures/input_fixture.tsv",
        "mags": False,
        "bins": True,
        "out": str(tmp_path),
        "force": False,
        "live": False,
        "test_suffix": "unit",
        "tpa": False,
        "centre_name": "EMG",
        "private": False,
    }

    uploader = GenomeUpload(args)

    xml_bytes = b"<SAMPLE_SET/>"

    single_path = uploader.write_genomes_xml(xml_bytes, 1, 1)
    batch_1_path = uploader.write_genomes_xml(xml_bytes, 1, 2)
    batch_2_path = uploader.write_genomes_xml(xml_bytes, 2, 2)
    retry_path = uploader.write_genomes_xml(xml_bytes, 2, 2, retry=True)

    assert single_path == uploader.samples_xml
    assert batch_1_path.name == "genome_samples_batch_1.xml"
    assert batch_2_path.name == "genome_samples_batch_2.xml"
    assert retry_path.name == "genome_samples_batch_2_retry.xml"
