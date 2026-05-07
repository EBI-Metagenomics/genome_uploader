import logging
import re
import time
import xml.dom.minidom as minidom

import requests
from retry import retry

from genomeuploader.ena import CredentialsManager
from genomeuploader.exceptions import EnaQueueTimeoutError

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

HTTP_ACCEPTED_CODE = 202
HTTP_TIMEOUT_CODES = (408, 504)

ENA_PROD_BASE_URL = "https://www.ebi.ac.uk/ena/submit/webin-v2"
ENA_DEV_BASE_URL = "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2"


def identify_registered_genomes(message):
    """
    Parses an ENA error message to check for already registered genomes.

    Args:
        message (str): Error message from ENA receipt containing alias and accession info.

    Returns:
        dict: Dictionary mapping genome alias to accession for already registered genomes.
    """
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
    def __init__(self, sample_xml, submission_receipt, number_of_genomes, live=False):
        self.sample_xml = sample_xml
        self.submission_receipt = submission_receipt
        self.live = live
        self.auth = CredentialsManager.get_credentials()
        self.number_of_genomes = number_of_genomes

    def poll_submission_receipt(self, poll_url: str, timeout_seconds: int = 600, poll_interval_seconds: int = 5) -> str:
        """
        Poll the ENA Webin queue endpoint until a final XML receipt is returned.

        The queue endpoint returns HTTP 202 while ENA is still processing the
        submission. In that case, the method waits for ``poll_interval_seconds``
        and retries until either a final receipt is available or the overall
        timeout is reached. HTTP 408 and 504 responses are treated as ENA queue
        timeouts and converted into ``EnaQueueTimeoutError``.

        Args:
            poll_url (str): Poll endpoint returned by the Webin submission API.
            timeout_seconds (int, optional): Maximum total time to wait for a
                final receipt before raising ``EnaQueueTimeoutError``.
                Defaults to 600 seconds (10 minutes).
            poll_interval_seconds (int, optional): Delay between polling
                attempts while the submission remains queued. Defaults to 5.

        Returns:
            str: Final XML receipt content returned by ENA.

        Raises:
            EnaQueueTimeoutError: If ENA returns a timeout status or if the
                final receipt is not available before the deadline.
        """
        deadline = time.monotonic() + timeout_seconds
        headers = {"Accept": "application/xml"}

        while time.monotonic() < deadline:
            poll_response = requests.get(poll_url, headers=headers, auth=self.auth)

            if poll_response.status_code == HTTP_ACCEPTED_CODE:
                logger.info("Submission is still being processed; waiting for final receipt...")
                time.sleep(poll_interval_seconds)
                continue

            if poll_response.status_code in HTTP_TIMEOUT_CODES:
                raise EnaQueueTimeoutError(
                    "ENA async submission timed out while processing. "
                    "This submission is likely too large and should be split into smaller batches. "
                    f"Polling response payload: {poll_response.text}"
                )

            poll_response.raise_for_status()
            return poll_response.text

        raise EnaQueueTimeoutError(
                    "ENA async submission timed out while processing. "
                    "This submission is likely too large and should be split into smaller batches. "
                    f"Polling response payload: {poll_response.text}"
                )

    def parse_receipt(self, receipt_content: str) -> dict:
        """
        Parse an ENA XML receipt and extract genome alias-to-accession mappings.

        Successful receipts return all submitted samples from the ``SAMPLE``
        elements. Failed receipts are inspected for ``ERROR`` elements that may
        still contain accessions for genomes already registered in ENA, and
        those previously registered genomes are returned when present.

        Args:
            receipt_content (str): XML receipt content returned by ENA.

        Returns:
            dict: Mapping of genome alias to accession for successfully parsed
                or previously registered genomes. Returns an empty dictionary if
                no accessions can be recovered from the receipt.
        """
        receipt_xml = minidom.parseString(receipt_content)
        receipt = receipt_xml.getElementsByTagName("RECEIPT")
        success = receipt[0].attributes["success"].value
        alias_dict = {}

        if success == "true":
            for sample in receipt_xml.getElementsByTagName("SAMPLE"):
                alias_dict[sample.attributes["alias"].value] = sample.attributes["accession"].value
            logger.info(f"{len(alias_dict)} genome samples successfully registered.")
            return alias_dict

        errors = receipt_xml.getElementsByTagName("ERROR")
        final_error = "".join(f"\n\t{error.firstChild.nodeValue.strip()}" for error in errors)

        registered_genomes = identify_registered_genomes(final_error)
        if registered_genomes:
            logger.info("Some previously submitted genomes were retrieved from the receipt")
            return registered_genomes

        logger.info("No previously submitted genomes retrieved from the receipt")
        return alias_dict

    @retry(
        exceptions=(requests.exceptions.RequestException),
        tries=3,
        delay=15,
        backoff=2,
        max_delay=120,
        logger=logger,
    )
    def register_genome_samples_in_ena(self):
        """
        Submits genome sample and submission XML files to ENA and parses the
        receipt for registration results.
        Handles both live and test modes, and extracts successfully registered
        genomes as well as previously registered ones from error messages.
        Updates and returns a dictionary mapping genome alias to accession for
        all registered genomes.
        
        The entire submission+polling workflow is retried with exponential backoff
        on transient errors (connection errors, HTTP 5xx errors).

        Returns:
            dict: Dictionary mapping genome alias to accession for all successfully or previously registered genomes.

        Raises:
            Exception: If the submission response does not contain expected fields.
        """
        base_url = ENA_PROD_BASE_URL if self.live else ENA_DEV_BASE_URL
        queue_url = f"{base_url}/submit/queue"
    
        mode = "live" if self.live else "test"
        logger.info(f"Registering genome samples using XML in {mode} mode.")

        submission_response = requests.post(
                queue_url,
                data=self.sample_xml.read_bytes(),
                headers={"Accept": "application/json", "Content-Type": "application/xml"},
                auth=self.auth,
            )
        submission_response.raise_for_status()

        queue_response = submission_response.json()
        submission_id = queue_response.get("submissionId")
        if not submission_id:
            raise Exception("ENA queue submission did not return a submissionId.")

        poll_url = queue_response.get("_links", {}).get("poll", {}).get("href")
        if not poll_url:
            raise Exception("ENA queue submission did not return a poll URL.")
        
        receipt_content = self.poll_submission_receipt(poll_url)

        # Write receipt XML to file for troubleshooting
        with open(self.submission_receipt, "w") as file:
            file.write(receipt_content)
            logger.info(f"Receipt XML written to {self.submission_receipt}")

        alias_dict = self.parse_receipt(receipt_content)
        if len(alias_dict) == self.number_of_genomes:
            logger.info("All genomes were registered")
        else:
            logger.info("For the re-registration some genomes will be excluded from the XML receipt.")
        return alias_dict
