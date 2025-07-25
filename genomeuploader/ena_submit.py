import logging
import xml.dom.minidom as minidom
import re

import requests

from genomeuploader.ena import CredentialsManager

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def identify_registered_genomes(message):
    alias_dict = {}
    pattern = r'alias: "([^"]+)"[^:]+accession: "([^"]+)"'
    for line in message.split("\n"):
        match = re.search(pattern, line)
        if match:
            alias = match.group(1)
            accession = match.group(2)
            alias_dict[alias] = accession
            logger.info(f"Found existing genome {alias} registered with {accession}")
    return alias_dict


class EnaSubmit:
    def __init__(self, sample_xml, submission_xml, number_of_genomes, live=False):
        self.sample_xml = sample_xml
        self.submission_xml = submission_xml
        self.live = live
        self.auth = CredentialsManager.get_credentials()
        self.number_of_genomes = number_of_genomes


    def handle_genomes_registration(self):
        live_sub, mode = "", "live"

        if not self.live:
            live_sub = "dev"
            mode = "test"

        url = f"https://www{live_sub}.ebi.ac.uk/ena/submit/drop-box/submit/"

        logger.info(f"Registering sample xml in {mode} mode.")

        try:
            with self.submission_xml.open("r") as sub_file, self.sample_xml.open("r") as sample_file:
                f = {"SUBMISSION": sub_file, "SAMPLE": sample_file}
                submission_response = requests.post(url, files=f, auth=self.auth)
                submission_response.raise_for_status()
        except requests.exceptions.RequestException as e:
            raise Exception(f"Failed to submit genomes to ENA: {e}")

        receipt_xml = minidom.parseString(submission_response.content.decode("utf-8"))
        receipt = receipt_xml.getElementsByTagName("RECEIPT")
        success = receipt[0].attributes["success"].value
        alias_dict = {}
        if success == "true":
            samples = receipt_xml.getElementsByTagName("SAMPLE")
            for s in samples:
                sra_acc = s.attributes["accession"].value
                alias = s.attributes["alias"].value
                alias_dict[alias] = sra_acc
            logger.info(f"{str(len(alias_dict))} genome samples successfully registered.")
        # check errors and search for existing accessions
        elif success == "false":
            errors = receipt_xml.getElementsByTagName("ERROR")
            final_error = ""
            for error in errors:
                final_error += f"\n\t{error.firstChild.nodeValue.strip()}"
            # check are there already registered genomes
            registered_genomes = identify_registered_genomes(final_error)
            if registered_genomes:
                alias_dict.update(registered_genomes)
                logger.info("Some previously submitted genomes were retrieved from the receipt")
            else:
                logger.info("No previously submitted genomes retrieved from the receipt")
        if len(alias_dict) == self.number_of_genomes:
            logger.info("All genomes were registered")
        else:
            logger.info("For the re-registration some genomes will be excluded from the XML receipt.")
        return alias_dict
