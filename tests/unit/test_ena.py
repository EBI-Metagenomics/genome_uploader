import pytest
import responses
from requests.exceptions import ConnectionError, HTTPError
import json

from genomeuploader.ena import EnaQuery
#from assembly_uploader.webin_utils import ENA_WEBIN, ENA_WEBIN_PASSWORD


def read_json(path):
    with open(path, 'r') as jpath:
        return json.load(jpath)

def read_xml(path):
    with open(path, 'r') as xpath:
        return xpath.read()


def test_ena_study(public_study_data, private_study_data, public_study_json, private_study_xml):
    responses.add(
        responses.POST,
        "https://www.ebi.ac.uk/ena/portal/api/search",
        json= read_json(public_study_json)
    )

    responses.add(
        responses.GET,
        "https://www.ebi.ac.uk/ena/submit/report/studies/xml/ERP125469",
        body= read_xml(private_study_xml),
        content_type='application/xml'
    )

    ena_study_public = EnaQuery(
        accession="ERP125469",
        query_type="study",
        private=False
    )

    ena_study_private = EnaQuery(
        accession="ERP125469",
        query_type="study",
        private=True,
    )
    assert ena_study_public.build_query() == public_study_data
    assert ena_study_private.build_query() == private_study_data


def test_ena_run(public_run_data, private_run_data, public_run_json, private_run_json):
    responses.add(
        responses.POST,
        "https://www.ebi.ac.uk/ena/portal/api/search",
        json = read_json(public_run_json)
    )

    responses.add(
        responses.GET,
        "https://www.ebi.ac.uk/ena/submit/report/runs/ERR4918394",
        json = read_json(private_run_json)
    )

    ena_run_public = EnaQuery(
        accession="ERR4918394",
        query_type="run",
        private=False
    )

    ena_run_private = EnaQuery(
        accession="ERR4918394",
        query_type="run",
        private=True
    )
    assert ena_run_public.build_query() == public_run_data
    assert ena_run_private.build_query() == private_run_data

def test_ena_run_from_assembly(public_run_from_assembly_xml, private_run_from_assembly_xml):
    responses.add(
        responses.GET,
        "https://www.ebi.ac.uk/ena/browser/api/xml/ERZ2626953",
        body = read_xml(public_run_from_assembly_xml),
        content_type='application/xml'
    )

    responses.add(
        responses.GET,
        "https://www.ebi.ac.uk/ena/submit/report/analyses/xml/ERZ2626953",
        body = read_xml(private_run_from_assembly_xml),
        content_type='application/xml'
    )

    ena_run_from_assembly_public = EnaQuery(
        accession="ERZ2626953",
        query_type="run_assembly",
        private=False
    )

    ena_run_from_assembly_public = EnaQuery(
        accession="ERZ2626953",
        query_type="run_assembly",
        private=True
    )

    assert ena_run_from_assembly_public.build_query() and ena_run_from_assembly_public.build_query() == "ERR4918394"

def test_ena_study_runs(public_study_runs_json, private_study_runs_json):
    """response mock limited to 10 runs for the study, need to check if there is page limit"""
    responses.add(
        responses.POST,
        "https://www.ebi.ac.uk/ena/portal/api/search",
        json = read_json(public_study_runs_json)
    )

    responses.add(
        responses.GET,
        "https://www.ebi.ac.uk/ena/submit/report/runs/ERP125469",
        json = read_json(private_study_runs_json)
    )

    ena_study_runs_public = EnaQuery(
        accession="ERP125469",
        query_type="study_runs",
        private=False
    )

    ena_study_runs_private = EnaQuery(
        accession="ERP125469",
        query_type="study_runs",
        private=True
    )

    assert len(ena_study_runs_public.build_query()) and len(ena_study_runs_private.build_query()) == 10 

def test_ena_sample(public_sample_data, private_sample_data, public_sample_json, private_sample_xml):
    responses.add(
        responses.POST,
        "https://www.ebi.ac.uk/ena/portal/api/search",
        json = read_json(public_sample_json)
    )

    responses.add(
        responses.GET,
        "https://www.ebi.ac.uk/ena/submit/report/samples/xml/ERS5444411",
        body = read_xml(private_sample_xml)
    )

    ena_sample_public = EnaQuery(
        accession="SAMEA7687881",
        query_type="sample",
        private=False
    )

    ena_sample_private = EnaQuery(
        accession="ERS5444411",
        query_type="sample",
        private=True
    )

    assert ena_sample_public.build_query() == public_sample_data
    assert ena_sample_private.build_query() == private_sample_data