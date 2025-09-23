#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2017-2025 EMBL - European Bioinformatics Institute
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

from genomeuploader.config import (
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
from genomeuploader.exceptions import (
    CredentialsNotSetError,
    EnaEmptyResponseError,
    EnaParseError,
    InvalidAccessionError,
)

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)


def get_default_connection_headers():
    """
    Returns default HTTP headers for ENA API requests.
    Returns:
        dict: Dictionary containing HTTP headers.
    """
    return {
        "headers": {
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "*/*",
        }
    }


def get_default_params():
    """
    Returns default parameters for ENA API requests.
    Returns:
        dict: Dictionary of default API parameters.
    """
    return {"format": "json", "includeMetagenomes": True, "dataPortal": "ena"}


def parse_accession(accession):
    """
    Determines the accession type based on its prefix.
    Args:
        accession (str): ENA accession string.
    Returns:
        str: Accession type (e.g., 'run', 'study', etc.).
    Raises:
        InvalidAccessionError: If the accession format is not recognized.
    """
    for prefix, acc_type in ACCESSION_MAP.items():
        if accession.startswith(prefix):
            return acc_type
    logger.error(f"Unrecognised accession format: '{accession}'")
    raise InvalidAccessionError("Invalid accession: {accession}")


def configure_credentials(env_filename=USER_ENV_FILE_PATH):
    """
    Loads ENA credentials from environment variables or .env files.
    Args:
        env_filename (str): Name of the environment file to search for.
    Returns:
        tuple: (username, password) for ENA authentication, or (None, None) if not found.
    """
    search_paths = [Path.home() / env_filename, Path.cwd() / env_filename, Path.cwd() / ".env"]

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
        return None, None

    return username, password


def query_taxid(taxid):
    """
    Queries ENA taxonomy API for a scientific name given a taxid.
    Args:
        taxid (str or int): NCBI taxonomic ID.
    Returns:
        str: Scientific name corresponding to the taxid, or empty string if not found.
    """
    url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{taxid}"
    response = requests.get(url)

    try:
        # Will raise exception if response status code is non-200
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        logger.error(f"Request failed {url} with error {e}")
        return False

    res = response.json()

    return res.get("scientificName", "")


def query_scientific_name(scientific_name, search_rank=False):
    """
    Queries ENA taxonomy API for a scientific name and optionally its rank.
    If the name exists, it checks if it's submittable and retrieves its taxid.
    Args:
        scientific_name (str): Scientific name to query.
        search_rank (bool): If True, also return taxonomic rank.
    Returns:
        tuple or bool: If search_rank is True, returns (submittable, taxid,
        rank). Otherwise, returns (submittable, taxid).
    """
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
            raise CredentialsNotSetError(
                "Credentials for ENA are not set. Please set ENA_WEBIN and ENA_WEBIN_PASSWORD environment variables."
            )
        return username, password


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
        """
        Sends a POST request to the ENA public API.
        Args:
            data (dict): Data payload for the POST request.
        Returns:
            requests.Response: Response object from the ENA API.
        """
        response = requests.post(self.public_url, data=data, **get_default_connection_headers())
        return response

    def get_request(self, url):
        """
        Sends a GET request to the ENA API (private or public).
        Args:
            url (str): URL to send the GET request to.
        Returns:
            requests.Response: Response object from the ENA API.
        """
        if self.private:
            response = requests.get(url, auth=self.auth)
        else:
            response = requests.get(url)
        return response

    def get_data_or_raise(self, response, mode="single_json"):
        """
        Parses the response from ENA and raises an error if the response
        is empty or invalid. It varies parsing based on the input format.
        Args:
            response (requests.Response): Response object from ENA API.
            mode (str): Parsing mode ('xml', 'full_json', 'single_json').
        Returns:
            object: Parsed response data (XML document or JSON object).
        Raises:
            EnaEmptyResponseError, EnaParseError: If response is empty or cannot be parsed.
        """
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
        """
        Retries a request to the ENA API up to RETRY_COUNT times on
        connection or timeout errors.
        Args:
            request (callable): Function to call for the request (get or post).
            *args: Arguments for the request function.
            **kwargs: Keyword arguments for the request function.
        Returns:
            requests.Response: Response object from the ENA API.
        Raises:
            ValueError, HTTPError, RequestException: If all retries fail or other errors occur.
        """
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
                logger.error(f"HTTP response has an error status: {http_err}")
                raise
            except RequestException as req_err:
                logger.error(f"Network-related error status: {req_err}")
                raise
            #   should hopefully encompass all other issues...
            except Exception as e:
                logger.error(f"An unexpected error occurred: {e}")
                raise

    def _fetch_ena_data(self, *, url=None, data=None, method="get", mode="single_json", reformatter=None):
        """
        Fetches and parses data from ENA using the specified method and mode.
        Args:
            url (str, optional): URL for GET requests.
            data (dict, optional): Data payload for POST requests.
            method (str): 'get' or 'post'.
            mode (str): Parsing mode ('xml', 'full_json', 'single_json').
            reformatter (callable, optional): Function to reformat parsed data.
        Returns:
            object: Reformatted or parsed response data from ENA.
        """
        request_func = self.get_request if method == "get" else self.post_request
        request_input = url or data
        response = self.retry_or_handle_request_error(request_func, request_input)
        try:
            parsed = self.get_data_or_raise(response, mode=mode)
        except (EnaEmptyResponseError, EnaParseError) as e:
            logger.error(e)
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
        logger.info(f"{self.accession} private run returned from ENA")
        return data

    def _get_public_run(self):
        data = get_default_params()
        data.update(
            {"result": "read_run", "query": f'run_accession="{self.accession}"', "fields": "secondary_study_accession,sample_accession"}
        )
        result = self._fetch_ena_data(data=data, method="post", mode="single_json")
        logger.info(f"{self.accession} public run returned from ENA")
        return result

    def _get_private_study(self):
        url = f"{self.private_url}/studies/xml/{self.accession}"

        def reformatter(xml_doc):
            desc = xml_doc.getElementsByTagName("STUDY_DESCRIPTION")[0].firstChild.nodeValue
            return {"study_accession": self.accession, "study_description": desc}

        result = self._fetch_ena_data(url=url, mode="xml", reformatter=reformatter)
        logger.info(f"{self.accession} private study returned from ENA")
        return result

    def _get_public_study(self):
        data = get_default_params()
        data.update({"result": "study", "query": f'{self.acc_type}="{self.accession}"', "fields": STUDY_DEFAULT_FIELDS})
        result = self._fetch_ena_data(data=data, method="post", mode="single_json")
        logger.info(f"{self.accession} public study returned from ENA")
        return result

    def _get_private_run_from_assembly(self):
        url = f"{self.private_url}/analyses/xml/{self.accession}"

        def reformatter(xml_doc):
            return xml_doc.getElementsByTagName("RUN_REF")[0].attributes["accession"].value

        result = self._fetch_ena_data(url=url, mode="xml", reformatter=reformatter)
        logger.info(f"private run from the assembly {self.accession} returned from ENA")
        return result

    def _get_public_run_from_assembly(self):
        url = f"{self.browser_url}/analyses/xml/{self.accession}"

        def reformatter(xml_doc):
            return xml_doc.getElementsByTagName("RUN_REF")[0].attributes["accession"].value

        result = self._fetch_ena_data(url=url, mode="xml", reformatter=reformatter)
        logger.info(f"public run from the assembly {self.accession} returned from ENA")
        return result

    def _get_private_study_runs(self):
        #   Pagination documentation unclear - offset not working. Using hardcoded max 2000 default. TO MODIFY
        url = f"{self.private_url}/runs/{self.accession}?format=json&max-results=2000"

        def reformatter(runs):
            result = []
            for run in runs:
                run_data = run["report"]
                updated_data = {RUN_KEY_RENAME_MAP.get(k, k): v for k, v in run_data.items()}
                result.append(updated_data)
            return result

        result = self._fetch_ena_data(url=url, mode="full_json", reformatter=reformatter)
        logger.info(f"found {len(result)} runs for study {self.accession}")
        logger.info(f"private runs from study {self.accession}, returned from ENA")
        return result

    def _get_public_study_runs(self):
        data = get_default_params()
        data.update({"result": "read_run", "fields": STUDY_RUN_DEFAULT_FIELDS, "query": f"{self.acc_type}={self.accession}"})
        result = self._fetch_ena_data(data=data, method="post", mode="full_json")
        logger.info(f"public runs from study {self.accession}, returned from ENA")
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
        logger.info(f"{self.accession} private sample returned from ENA")
        return result

    def _get_public_sample(self):
        data = get_default_params()
        data.update({"result": "sample", "fields": SAMPLE_DEFAULT_FIELDS, "query": f"{self.acc_type}={self.accession}"})
        result = self._fetch_ena_data(data=data, method="post", mode="single_json")
        logger.info(f"{self.accession} public sample returned from ENA")
        return result

    def build_query(self):
        """
        Executes the appropriate ENA query (private or public) for the specified
        query type. If the private flag is given, assume private data and try
        private APIs. ENA also has cases where a run may be private but the
        sample may be public etc. Hence always try public if private fails
        Returns:
            dict or list: ENA response data for the requested query type (study, run, sample, etc.).
        """
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
                logger.info(f"Private API for {self.query_type} failed, trying public.")
                ena_response = public_api()
        else:
            ena_response = public_api()

        return ena_response
