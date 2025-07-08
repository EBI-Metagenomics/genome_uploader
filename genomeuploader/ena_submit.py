import logging
import requests
import xml.dom.minidom as minidom


from genomeuploader.ena import configure_credentials

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class EnaSubmit:
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

        logger.info("Registering sample xml in {} mode.".format(mode))

        try:
            with open(self.submission_xml, "r") as sub_file, open(self.sample_xml, "r") as sample_file:
                f = {"SUBMISSION": sub_file, "SAMPLE": sample_file}
                submission_response = requests.post(url, files=f, auth=self.auth)
                submission_response.raise_for_status()
        except requests.exceptions.RequestException as e:
            raise Exception(f"Failed to submit genomes to ENA: {e}")

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

        logger.info("{} genome samples successfully registered.".format(str(len(alias_dict))))

        return alias_dict