#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2017-2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import json
import logging
import os
import xml.dom.minidom as minidom
from pathlib import Path
from time import sleep

import requests
from dotenv import load_dotenv
from requests.exceptions import ConnectionError, HTTPError, RequestException, Timeout

from .config import (
    ACCESSION_MAP,
    ENA_BROWSER_URL,
    ENA_SEARCH_URL,
    ENA_SUBMISSION_URL,
    RETRY_COUNT,
    RUN_KEY_RENAME_MAP,
    SAMPLE_ATTRIBUTE_MAP,
    SAMPLE_DEFAULT_FIELDS,
    STUDY_DEFAULT_FIELDS,
    STUDY_RUN_DEFAULT_FIELDS,
    USER_ENV_FILE_PATH,
)
from .exceptions import EnaEmptyResponseError, EnaParseError, InvalidAccessionError

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)


def get_default_connection_headers():
    return {
        "headers": {
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "*/*",
        }
    }


def get_default_params():
    return {"format": "json", "includeMetagenomes": True, "dataPortal": "ena"}


def parse_accession(accession):
    for prefix, acc_type in ACCESSION_MAP.items():
        if accession.startswith(prefix):
            return acc_type
    logging.error(f"Unrecognized accession format: '{accession}'")
    raise InvalidAccessionError("Invalid accession: {accession}")


def configure_credentials(env_filename=USER_ENV_FILE_PATH):
    search_paths = [Path.home() / env_filename, Path.cwd() / env_filename, Path.cwd / ".env"]

    for env_path in search_paths:
        if env_path.exists():
            logger.debug(f"Loading environment variables from {env_path}")
            load_dotenv(dotenv_path=env_path)
            break
    else:
        logger.debug("No environment file found; relying on system environment variables.")

    username = os.getenv("ENA_WEBIN")
    password = os.getenv("ENA_WEBIN_PASSWORD")

    if not username or not password:
        logger.warning("ENA_WEBIN and/or ENA_WEBIN_PASSWORD not found in environment variables.")

    return username, password


def query_taxid(taxid):
    url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{taxid}"
    response = requests.get(url)

    try:
        # Will raise exception if response status code is non-200
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        logging.error(f"Request failed {url} with error {e}")
        return False

    res = response.json()

    return res.get("scientificName", "")


def query_scientific_name(scientific_name, search_rank=False):
    url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{scientific_name}"
    response = requests.get(url)

    try:
        # Will raise exception if response status code is non-200
        response.raise_for_status()
    except requests.exceptions.HTTPError:
        if search_rank:
            return False, "", ""
        else:
            return False, ""

    try:
        res = response.json()[0]
    except IndexError:
        if search_rank:
            return False, "", ""
        else:
            return False, ""

    submittable = res.get("submittable", "").lower() == "true"
    taxid = res.get("taxId", "")
    rank = res.get("rank", "")

    if search_rank:
        return submittable, taxid, rank
    else:
        return submittable, taxid


class CredentialsManager:
    @staticmethod
    def get_credentials():
        username, password = configure_credentials()
        if not username or not password:
            logging.error("ENA_WEBIN and ENA_WEBIN_PASSWORD are not set.")
            return None
        return (username, password)


class EnaQuery:
    def __init__(self, accession, query_type, private=False):
        self.private_url = ENA_SUBMISSION_URL
        self.public_url = ENA_SEARCH_URL
        self.browser_url = ENA_BROWSER_URL
        self.accession = accession
        self.acc_type = parse_accession(accession)
        self.auth = CredentialsManager.get_credentials()
        self.private = private
        self.query_type = query_type

    def post_request(self, data):
        response = requests.post(self.public_url, data=data, **get_default_connection_headers())
        return response

    def get_request(self, url):
        if self.private:
            response = requests.get(url, auth=self.auth)
        else:
            response = requests.get(url)
        return response

    def get_data_or_raise(self, response, mode="single_json"):
        data_txt = response.text.strip()
        if not data_txt:
            raise EnaEmptyResponseError(f"{self.accession}: Empty or missing response text.")

        parsers = {
            "xml": lambda txt: minidom.parseString(txt),
            "full_json": lambda txt: json.loads(txt),
            "single_json": lambda txt: json.loads(txt)[0],
        }

        if mode not in parsers:
            raise EnaParseError(f"Unknown mode: {mode}")

        try:
            return parsers[mode](data_txt)
        except Exception as e:
            raise EnaParseError(f"{self.accession}: Failed to parse response as {mode}: {e}")

    def retry_or_handle_request_error(self, request, *args, **kwargs):
        attempt = 0
        while attempt < RETRY_COUNT:
            try:
                response = request(*args, **kwargs)
                response.raise_for_status()
                return response
            #   all other RequestExceptions are raised below
            except (Timeout, ConnectionError) as retry_err:
                attempt += 1
                if attempt >= RETRY_COUNT:
                    raise ValueError(f"Could not find {self.accession} in ENA after {attempt} attempts. Error: {retry_err}")
                sleep(1)
            except HTTPError as http_err:
                logging.error(f"HTTP response has an error status: {http_err}")
                raise
            except RequestException as req_err:
                logging.error(f"Network-related error status: {req_err}")
                raise
            #   should hopefully encompass all other issues...
            except Exception as e:
                logging.error(f"An unexpected error occurred: {e}")
                raise

    def _fetch_ena_data(self, *, url=None, data=None, method="get", mode="single_json", reformatter=None):
        request_func = self.get_request if method == "get" else self.post_request
        request_input = url or data
        response = self.retry_or_handle_request_error(request_func, request_input)
        try:
            parsed = self.get_data_or_raise(response, mode=mode)
        except (EnaEmptyResponseError, EnaParseError) as e:
            logging.error(e)
            return None
        return reformatter(parsed) if reformatter else parsed

    def _get_private_run(self):
        url = f"{self.private_url}/runs/{self.accession}"

        def reformatter(run):
            run_data = run["report"]
            return {
                "run_accession": run_data["id"],
                "secondary_study_accession": run_data["studyId"],
                "sample_accession": run_data["sampleId"],
            }

        data = self._fetch_ena_data(url=url, mode="single_json", reformatter=reformatter)
        logging.info(f"{self.accession} private run returned from ENA")
        return data

    def _get_public_run(self):
        data = get_default_params()
        data.update(
            {"result": "read_run", "query": f'run_accession="{self.accession}"', "fields": "secondary_study_accession,sample_accession"}
        )
        result = self._fetch_ena_data(data=data, method="post", mode="single_json")
        logging.info(f"{self.accession} public run returned from ENA")
        return result

    def _get_private_study(self):
        url = f"{self.private_url}/studies/xml/{self.accession}"

        def reformatter(xml_doc):
            desc = xml_doc.getElementsByTagName("STUDY_DESCRIPTION")[0].firstChild.nodeValue
            return {"study_accession": self.accession, "study_description": desc}

        result = self._fetch_ena_data(url=url, mode="xml", reformatter=reformatter)
        logging.info(f"{self.accession} private study returned from ENA")
        return result

    def _get_public_study(self):
        data = get_default_params()
        data.update({"result": "study", "query": f'{self.acc_type}="{self.accession}"', "fields": STUDY_DEFAULT_FIELDS})
        result = self._fetch_ena_data(data=data, method="post", mode="single_json")
        logging.info(f"{self.accession} public study returned from ENA")
        return result

    def _get_private_run_from_assembly(self):
        url = f"{self.private_url}/analyses/xml/{self.accession}"

        def reformatter(xml_doc):
            return xml_doc.getElementsByTagName("RUN_REF")[0].attributes["accession"].value

        result = self._fetch_ena_data(url=url, mode="xml", reformatter=reformatter)
        logging.info(f"private run from the assembly {self.accession} returned from ENA")
        return result

    def _get_public_run_from_assembly(self):
        url = f"{self.browser_url}/analyses/xml/{self.accession}"

        def reformatter(xml_doc):
            return xml_doc.getElementsByTagName("RUN_REF")[0].attributes["accession"].value

        result = self._fetch_ena_data(url=url, mode="xml", reformatter=reformatter)
        logging.info(f"public run from the assembly {self.accession} returned from ENA")
        return result

    def _get_private_study_runs(self):
        #   Pagination documentation unclear - offest not working. Using hardcoded max 2000 default. TO MODIFY
        url = f"{self.private_url}/runs/{self.accession}?format=json&max-results=2000"

        def reformatter(runs):
            result = []
            for run in runs:
                run_data = run["report"]
                updated_data = {RUN_KEY_RENAME_MAP.get(k, k): v for k, v in run_data.items()}
                result.append(updated_data)
            return result

        result = self._fetch_ena_data(url=url, mode="full_json", reformatter=reformatter)
        logging.info(f"found {len(result)} runs for study {self.accession}")
        logging.info(f"private runs from study {self.accession}, returned from ENA")
        return result

    def _get_public_study_runs(self):
        data = get_default_params()
        data.update({"result": "read_run", "fields": STUDY_RUN_DEFAULT_FIELDS, "query": f"{self.acc_type}={self.accession}"})
        result = self._fetch_ena_data(data=data, method="post", mode="full_json")
        logging.info(f"public runs from study {self.accession}, returned from ENA")
        return result

    def _get_private_sample(self):
        url = f"{self.private_url}/samples/xml/{self.accession}"

        def reformatter(xml_doc):
            reformatted = {"sample_accession": self.accession}
            for attr in xml_doc.getElementsByTagName("SAMPLE_ATTRIBUTE"):
                tag = attr.getElementsByTagName("TAG")[0].firstChild.nodeValue
                val = attr.getElementsByTagName("VALUE")[0].firstChild.nodeValue
                field = SAMPLE_ATTRIBUTE_MAP.get(tag)
                if field:
                    reformatted[field] = val
            return reformatted

        result = self._fetch_ena_data(url=url, mode="xml", reformatter=reformatter)
        logging.info(f"{self.accession} private sample returned from ENA")
        return result

    def _get_public_sample(self):
        data = get_default_params()
        data.update({"result": "sample", "fields": SAMPLE_DEFAULT_FIELDS, "query": f"{self.acc_type}={self.accession}"})
        result = self._fetch_ena_data(data=data, method="post", mode="single_json")
        logging.info(f"{self.accession} public sample returned from ENA")
        return result

    def build_query(self):
        """If the private flag is given, assume private data and try private APIs.
        ENA also has cases where a run may be private but the sample may be public etc. Hence always try
        public if private fails"""
        api_map = {
            "study": (self._get_private_study, self._get_public_study),
            "run": (self._get_private_run, self._get_public_run),
            "run_assembly": (self._get_private_run_from_assembly, self._get_public_run_from_assembly),
            "study_runs": (self._get_private_study_runs, self._get_public_study_runs),
            "sample": (self._get_private_sample, self._get_public_sample),
        }

        private_api, public_api = api_map.get(self.query_type, (None, None))

        if self.private:
            try:
                ena_response = private_api()
            except Exception:
                logging.info(f"Private API for {self.query_type} failed, trying public.")
                ena_response = public_api()
        else:
            ena_response = public_api()

        return ena_response
