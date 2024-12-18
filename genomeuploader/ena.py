
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


import requests
import json
import logging
import os
import sys 
from time import sleep
from pathlib import Path
from dotenv import load_dotenv

import xml.dom.minidom as minidom
from requests.exceptions import ConnectionError, HTTPError, RequestException, Timeout

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)


class NoDataException(ValueError):
    pass


RUN_DEFAULT_FIELDS = ','.join([
    'secondary_study_accession',
    'sample_accession'
])

STUDY_RUN_DEFAULT_FIELDS = ','.join([
    'sample_accession',
    'run_accession',
    'instrument_model'
])

ASSEMBLY_DEFAULT_FIELDS = 'sample_accession'


SAMPLE_DEFAULT_FIELDS = ','.join([
    'sample_accession',
    'secondary_sample_accession',
    'collection_date',
    'country',
    'location'
])

STUDY_DEFAULT_FIELDS = 'study_description'

RETRY_COUNT = 3

def get_default_connection_headers():
    return {
        "headers": {
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "*/*",
        }
    }

def get_default_params():
    return {
        'format': 'json',
        'includeMetagenomes': True,
        'dataPortal': 'ena'
    }

def parse_accession(accession):
    if accession.startswith("PRJ"):
        return "study_accession"
    elif "RP" in accession:
        return "secondary_study_accession"
    elif "RR" in accession:
        return "run_accession"
    elif "SAM" in accession:
        return "secondary_sample_accession"
    elif "RS" in accession:
        return "sample_accession"
    else:
        logging.error(f"{accession} is not a valid accession")
        sys.exit()

def configure_credentials():
    # Config file
    user_config = Path.home() / ".genome_uploader.config.env"
    if user_config.exists():
        logger.debug("Loading the env variables from ".format(user_config))
        load_dotenv(str(user_config))
    else:
        cwd_config = Path.cwd() / ".genome_uploader.config.env"
        if cwd_config.exists():
            logger.debug("Loading the variables from the current directory.")
            load_dotenv(str(cwd_config))
        else:
            logger.debug("Trying to load env variables from the .env file")
            # from a local .env file
            load_dotenv()

    username = os.getenv("ENA_WEBIN")
    password = os.getenv("ENA_WEBIN_PASSWORD")

    return username, password

def query_taxid(taxid):
    url = "https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{}".format(taxid)
    response = requests.get(url)

    try:
        # Will raise exception if response status code is non-200 
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        print("Request failed {} with error {}".format(url, e))
        return False
    
    res = response.json()
    
    return res.get("scientificName", "")

def query_scientific_name(scientific_name, search_rank=False):
    url = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{}".format(scientific_name)
    response = requests.get(url)
    
    try:
        # Will raise exception if response status code is non-200 
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
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



class EnaQuery():
    def __init__(self, accession, query_type, private=False):
        self.private_url = "https://www.ebi.ac.uk/ena/submit/report"
        self.public_url = "https://www.ebi.ac.uk/ena/portal/api/search"
        self.browser_url = "https://www.ebi.ac.uk/ena/browser/api/xml"
        self.accession = accession
        self.acc_type = parse_accession(accession)
        username, password = configure_credentials()
        if username is None or password is None:
            logging.error("ENA_WEBIN and ENA_WEBIN_PASSWORD are not set")
        if username and password:
            self.auth = (username, password)
        else:
            self.auth = None
        self.private = private
        self.query_type = query_type

    def post_request(self, data):
        response = requests.post(
            self.public_url, data=data, **get_default_connection_headers()
        )
        return response

    def get_request(self, url):
        if self.private:
            response = requests.get(url, auth=self.auth)
        else:
            response = requests.get(url)
        return response

    def get_data_or_handle_error(self, response, xml=False, full_json=False):
        try:
            data_txt = response.text
            if data_txt is None:
                if self.private:
                    logging.error(
                        f"{self.accession} private data is not present in the specified Webin account"
                    )
                else:
                    logging.error(f"{self.accession} public data does not exist")
            else:
                if xml:
                    return  minidom.parseString(data_txt)
                elif full_json: 
                    #   return the full json when more than one entry expected
                    return response.json()
                else:
                    #   only return the fist element in the list is the default. In these cases there should only be one entry returned 
                    return json.loads(response.txt)[0]
        except (IndexError, TypeError, ValueError, KeyError):
            logging.error(
                f"Failed to fetch {self.accession}, returned error: {response.text}"
            )     

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
                    raise ValueError(
                        f"Could not find {self.accession} in ENA after {RETRY_COUNT} attempts. Error: {retry_err}"
                    )
                sleep(1)
            except HTTPError as http_err:
                print(f"HTTP response has an error status: {http_err}")
                raise
            except RequestException as req_err:
                print(f"Network-related error status: {req_err}")
                raise
            #   should hopefully encompass all other issues...
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                raise

    def _get_private_run(self):
        url = f"{self.private_url}/runs/{self.accession}"
        response = self.retry_or_handle_request_error(self.get_request, url)
        run = self.get_data_or_handle_error(response)
        run_data = run["report"]
        reformatted_data = {
            "secondary_study_accession": run_data['studyID'],
            "sample_accession": run_data["sampleId"],
        }
        logging.info(f"{self.accession} private run returned from ENA")
        return reformatted_data

    def _get_public_run(self):
        data = get_default_params()
        data.update(
            {
            "result": "read_run",
            "query": f'run_accession="{self.accession}"',
            "fields": "secondary_study_accession,sample_accession",
        }
        )
        response = self.retry_or_handle_request_error(self.post_request, data)
        run = self.get_data_or_handle_error(response)
        logging.info(f"{self.accession} public run returned from ENA")
        return run   

    
    def _get_private_study(self):
        url = f"{self.private_url}/studies/xml/{self.accession}"
        response = self.retry_or_handle_request_error(self.get_request, url)
        study = self.get_data_or_handle_error(response, xml=True)
        study_desc = study.getElementsByTagName("STUDY_DESCRIPTION")[0].firstChild.nodeValue
        reformatted_data = {
            "study_description" : study_desc
        }
        logging.info(f"{self.accession} private study returned from ENA")
        return reformatted_data

    def _get_public_study(self):
        data = get_default_params()
        data.update(
            {
            "result": "study",
            "query": f'{self.acc_type}="{self.accession}"',
            "fields": STUDY_DEFAULT_FIELDS,
        }
        )
        response = self.retry_or_handle_request_error(self.post_request, data)
        study = self.get_data_or_handle_error(response)
        logging.info(f"{self.accession} public study returned from ENA")
        return study   

    def _get_private_run_from_assembly(self):
        url = f"{self.private_url}/analyses/xml/{self.accession}"
        reponse = self.retry_or_handle_request_error(self.get_request, url)
        parsed_xml = self.get_data_or_handle_error(reponse, xml=True)
        run_ref = parsed_xml.getElementsByTagName("RUN_REF")
        run = run_ref[0].attributes["accession"].value
        logging.info(f"private run from the assembly {self.accession} returned from ENA")
        return run

    def _get_public_run_from_assembly(self):
        url = f"{self.browser_url}/{self.accession}"
        reponse = self.retry_or_handle_request_error(self.get_request, url)
        parsed_xml = self.get_data_or_handle_error(reponse, xml=True)
        run_ref = parsed_xml.getElementsByTagName("RUN_REF")
        run = run_ref[0].attributes["accession"].value
        logging.info(f"public run from the assembly {self.accession} returned from ENA")

        return run  

    def _get_private_study_runs(self):
        url = f"{self.private_url}/runs/{self.accession}"
        response = self.retry_or_handle_request_error(self.get_request, url)
        runs = self.get_data_or_handle_error(response, full_json=True)
        reformatted_data = []
        for run in runs:
            run_data = run['report']
            if 'sampleId' in run_data:
                run_data['sample_accession'] = run_data.pop('sampleId')
            if 'id' in run_data:
                run_data['run_accession'] = run_data.pop('id')
            if 'instrumentModel' in run_data:
                run_data['instrument_model'] = run_data.pop('instrumentModel')
            reformatted_data.append(run_data)
        logging.info(f"private runs from study {self.accession}, returned from ENA")
        return reformatted_data


    def _get_public_study_runs(self):
        data = get_default_params()
        data.update(
            {
                "result": "read_run",
                "fields": STUDY_RUN_DEFAULT_FIELDS,
                "query": f"{self.acc_type}={self.accession}"
            }
        )
        response = self.retry_or_handle_request_error(self.post_request, data)
        runs = self.get_data_or_handle_error(response, full_json=True)
        logging.info(f"public runs from study {self.accession}, returned from ENA")
        return runs
    
    def _get_private_sample(self):
        url = f"{self.private_url}/samples/xml/{self.accession}"
        reponse = self.retry_or_handle_request_error(self.get_request, url)
        reformatted_data = {}
        parsed_xml = self.get_data_or_handle_error(reponse, xml=True)
        sample_attributes = parsed_xml.getElementsByTagName('SAMPLE_ATTRIBUTE')
        for attribute in sample_attributes:
            tag = attribute.getElementsByTagName('TAG')[0].firstChild.nodeValue
            if tag == "geographic location (country and/or sea)":
                reformatted_data['country'] = attribute.getElementsByTagName('VALUE')[0].firstChild.nodeValue
            if tag == "geographic location (latitude)":
                reformatted_data['latitude'] = attribute.getElementsByTagName('VALUE')[0].firstChild.nodeValue
            if tag == "geographic location (longitude)":
                reformatted_data['longitude'] = attribute.getElementsByTagName('VALUE')[0].firstChild.nodeValue
            if tag == "collection date":
                reformatted_data['collection_date'] = attribute.getElementsByTagName('VALUE')[0].firstChild.nodeValue
        logging.info(f"{self.accession} private sample returned from ENA")        
        return reformatted_data


    def _get_public_sample(self):
        data = get_default_params()
        data.update(
            {
                "result": "sample",
                "fields": SAMPLE_DEFAULT_FIELDS,
                "query":  f"{self.acc_type}={self.accession}"
            }
        )
        response = self.retry_or_handle_request_error(self.post_request, data)
        sample = self.get_data_or_handle_error(response)
        logging.info(f"{self.accession} public sample returned from ENA")
        return sample         


    def build_query(self):
        if self.query_type == "study":
            if self.private:
                ena_response = self._get_private_study()
            else:
                ena_response = self._get_public_study()
        elif self.query_type == "run":
            if self.private:
                ena_response = self._get_private_run()
            else:
                ena_response = self._get_public_run()
        elif self.query_type == "run_assembly":
            if self.private:
                ena_response = self._get_private_run_from_assembly()
            else:
                ena_response = self._get_public_run_from_assembly()
        elif self.query_type == "study_runs":
            if self.private:
                ena_response = self._get_private_study_runs()
            else:
                ena_response = self._get_public_study_runs()
        elif self.query_type == "sample":
            if self.private:
                ena_response = self._get_private_sample()
            else:
                ena_response = self._get_public_sample()
        return ena_response


class EnaSubmit():
    def __init__(self, sample_xml, submission_xml, live=False):
        self.sample_xml = sample_xml
        self.submission_xml = submission_xml
        self.live = live
        username, password = configure_credentials()
        if username is None or password is None:
            logging.error("ENA_WEBIN and ENA_WEBIN_PASSWORD are not set")
        if username and password:
            self.auth = (username, password)
        else:
            self.auth = None
        

    def handle_genomes_registration(self):
        live_sub, mode = "", "live"

        if not self.live:
            live_sub = "dev"
            mode = "test"

        url = "https://www{}.ebi.ac.uk/ena/submit/drop-box/submit/".format(live_sub)

        logger.info('Registering sample xml in {} mode.'.format(mode))

        f = {
            'SUBMISSION': open(self.submission_xml, 'r'),
            'SAMPLE': open(self.sample_xml, 'r')
        }

        submission_response = requests.post(url, files = f, auth = self.auth)

        if submission_response.status_code != 200:
            if str(submission_response.status_code).startswith('5'):
                raise Exception("Genomes could not be submitted to ENA as the server " +
                    "does not respond. Please again try later.")
            else:
                raise Exception("Genomes could not be submitted to ENA. HTTP response: " +
                    submission_response.reason)

        receipt_xml = minidom.parseString((submission_response.content).decode("utf-8"))
        receipt = receipt_xml.getElementsByTagName("RECEIPT")
        success = receipt[0].attributes["success"].value
        if success == "true":
            alias_dict = {}
            samples = receipt_xml.getElementsByTagName("SAMPLE")
            for s in samples:
                sra_acc = s.attributes["accession"].value
                alias = s.attributes["alias"].value
                alias_dict[alias] = sra_acc
        elif success == "false":
            errors = receipt_xml.getElementsByTagName("ERROR")
            final_error = "\tSome genomes could not be submitted to ENA. Please, check the errors below."
            for error in errors:
                final_error += "\n\t" + error.firstChild.data
            final_error += "\n\tIf you wish to validate again your data and metadata, "
            final_error += "please use the --force option."
            raise Exception(final_error)
        
        logger.info('{} genome samples successfully registered.'.format(str(len(alias_dict))))

        return alias_dict
    


        
