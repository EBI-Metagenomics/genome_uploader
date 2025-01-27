
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
from time import sleep

import xml.dom.minidom as minidom

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)


class NoDataException(ValueError):
    pass


RUN_DEFAULT_FIELDS = ','.join([
    'study_accession',
    'secondary_study_accession',
    'instrument_model',
    'run_accession',
    'sample_accession'
])

ASSEMBLY_DEFAULT_FIELDS = 'sample_accession'

SAMPLE_DEFAULT_FIELDS = ','.join([
    'sample_accession',
    'secondary_sample_accession',
    'collection_date',
    'country',
    'location'
])

STUDY_DEFAULT_FIELDS = ','.join([
    'study_accession',
    'secondary_study_accession',
    'study_description',
    'study_title'
])

RETRY_COUNT = 5


class ENA():
    def get_default_params(self):
        return {
            'format': 'json',
            'includeMetagenomes': True,
            'dataPortal': 'ena'
        }

    def post_request(self, data, webin, password):
        url = "https://www.ebi.ac.uk/ena/portal/api/search"
        auth = (webin, password)
        default_connection_headers = {
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "*/*"
        }
        response = requests.post(url, data=data, auth=auth, headers=default_connection_headers)
        
        return response

    def get_run(self, run_accession, webin, password, attempt=0, search_params=None):
        data = self.get_default_params()
        data['result'] = 'read_run'
        data['fields'] = RUN_DEFAULT_FIELDS
        data['query'] = 'run_accession=\"{}\"'.format(run_accession)

        if search_params:
            data.update(search_params)

        response = self.post_request(data, webin, password)
        
        if not response.ok and attempt > 2:
            raise ValueError("Could not retrieve run with accession {}, returned "
                "message: {}".format(run_accession, response.text))
        elif response.status_code == 204:
            if attempt < 2:
                attempt += 1
                sleep(1)
                return self.get_run(run_accession, webin, password, attempt)
            else:
                raise ValueError("Could not find run {} in ENA after {}"
                    " attempts".format(run_accession, RETRY_COUNT))
        try:
            run = json.loads(response.text)[0]
        except (IndexError, TypeError, ValueError):
            raise ValueError("Could not find run {} in ENA.".format(run_accession))
        except:
            raise Exception("Could not query ENA API for run {}: {}".format(run_accession, response.text))

        return run

    def get_run_from_assembly(self, assembly_name):
        manifestXml = minidom.parseString(requests.get("https://www.ebi.ac.uk" +
            "/ena/browser/api/xml/" + assembly_name).text)

        run_ref = manifestXml.getElementsByTagName("RUN_REF")
        run = run_ref[0].attributes["accession"].value
        
        return run

    def get_study(self, webin, password, study_accession):
        data = self.get_default_params()
        data['result'] = "study"
        data['fields'] = STUDY_DEFAULT_FIELDS
        data['query'] = 'study_accession="{}" OR secondary_study_accession="{}"' \
            .format(study_accession, study_accession)

        data['dataPortal'] = "ena"

        try:
            response = self.post_request(data, webin, password)
            if response.status_code == 204:
                raise NoDataException()
            try:
                studyList = response.json()
                assert len(studyList) == 1
                study = studyList[0]
            except (IndexError, TypeError, ValueError, KeyError) as e:
                raise e
            return study
        except NoDataException:
            print("No info found to fetch study {}".format(study_accession))
        except (IndexError, TypeError, ValueError, KeyError):
            print("Failed to fetch study {}, returned error: {}".format(study_accession, response.text))

        raise ValueError('Could not find study {} in ENA.'.format(study_accession))

    def get_study_runs(self, study_acc, webin, password, fields=None, search_params=None):
        data = self.get_default_params()
        data['result'] = 'read_run'
        data['fields'] = fields or RUN_DEFAULT_FIELDS
        data['query'] = '(study_accession=\"{}\" OR secondary_study_accession=\"{}\")'.format(study_acc, study_acc)

        if search_params:
            data.update(search_params)

        response = self.post_request(data, webin, password)
        
        if not response.ok:
            raise ValueError("Could not retrieve runs for study %s.", study_acc)
        
        if response.status_code == 204:
            return []

        try:
            runs = response.json()
        except:
            raise ValueError("Query against ENA API did not work. Returned "
                "message: {}".format(response.text))

        return runs

    def get_sample(self, sample_accession, webin, password, fields=None, search_params=None, attempt=0):
        data = self.get_default_params()
        data['result'] = 'sample'
        data['fields'] = fields or SAMPLE_DEFAULT_FIELDS
        data['query'] = ('(sample_accession=\"{acc}\" OR secondary_sample_accession'
            '=\"{acc}\") ').format(acc=sample_accession)

        if search_params:
            data.update(search_params)

        response = self.post_request(data, webin, password)
        
        if response.status_code == 200:
            sample = response.json()
            assert len(sample) == 1 
            return sample[0]

        if response.status_code == 204:
            if attempt < 2:
                new_params = {'dataPortal': 'metagenome' if data['dataPortal'] == 'ena' else 'ena'}
                attempt += 1
                return self.get_sample(sample_accession, webin, password, fields=fields, 
                    search_params=new_params, attempt=attempt)
            else:
                raise ValueError("Could not find sample {} in ENA after "
                    "{} attempts.".format(sample_accession, RETRY_COUNT))
        else:
            raise ValueError("Could not retrieve sample with accession {}. "
                "Returned message: {}".format(sample_accession, response.text))

    def query_taxid(self, taxid):
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

    def query_scientific_name(self, scientificName, searchRank=False):
        url = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{}".format(scientificName)
        response = requests.get(url)
        
        try:
            # Will raise exception if response status code is non-200 
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            if searchRank:
                return False, "", ""
            else:
                return False, ""
        
        try:
            res = response.json()[0]
        except IndexError:
            if searchRank:
                return False, "", ""
            else:
                return False, ""

        submittable = res.get("submittable", "").lower() == "true"
        taxid = res.get("taxId", "")
        rank = res.get("rank", "")

        if searchRank:
            return submittable, taxid, rank
        else:
            return submittable, taxid

    def handle_genomes_registration(self, sample_xml, submission_xml, webin, password, live=False):
        liveSub, mode = "", "live"

        if not live:
            liveSub = "dev"
            mode = "test"

        url = "https://www{}.ebi.ac.uk/ena/submit/drop-box/submit/".format(liveSub)

        logger.info('Registering sample xml in {} mode.'.format(mode))

        f = {
            'SUBMISSION': open(submission_xml, 'r'),
            'SAMPLE': open(sample_xml, 'r')
        }

        submissionResponse = requests.post(url, files = f, auth = (webin, password))

        if submissionResponse.status_code != 200:
            if str(submissionResponse.status_code).startswith('5'):
                raise Exception("Genomes could not be submitted to ENA as the server " +
                    "does not respond. Please again try later.")
            else:
                raise Exception("Genomes could not be submitted to ENA. HTTP response: " +
                    submissionResponse.reason)

        receiptXml = minidom.parseString((submissionResponse.content).decode("utf-8"))
        receipt = receiptXml.getElementsByTagName("RECEIPT")
        success = receipt[0].attributes["success"].value
        if success == "true":
            aliasDict = {}
            samples = receiptXml.getElementsByTagName("SAMPLE")
            for s in samples:
                sraAcc = s.attributes["accession"].value
                alias = s.attributes["alias"].value
                aliasDict[alias] = sraAcc
        elif success == "false":
            errors = receiptXml.getElementsByTagName("ERROR")
            finalError = "\tSome genomes could not be submitted to ENA. Please, check the errors below."
            for error in errors:
                finalError += "\n\t" + error.firstChild.data
            finalError += "\n\tIf you wish to validate again your data and metadata, "
            finalError += "please use the --force option."
            raise Exception(finalError)
        
        logger.info('{} genome samples successfully registered.'.format(str(len(aliasDict))))

        return aliasDict
