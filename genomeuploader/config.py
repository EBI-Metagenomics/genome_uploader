ACCESSION_MAP = {
    "PRJ": "study_accession",
    "ERP": "secondary_study_accession",
    "DRP": "secondary_study_accession",
    "SRP": "secondary_study_accession",
    "ERR": "run_accession",
    "DRR": "run_accession",
    "SRR": "run_accession",
    "SAM": "sample_accession",
    "ERS": "secondary_sample_accession",
    "DRS": "secondary_sample_accession",
    "SRS": "secondary_sample_accession",
    "ERZ": "analysis_accession",
    "DRZ": "analysis_accession",
    "SRZ": "analysis_accession",
}
USER_ENV_FILE_PATH = ".genome_uploader.config.env"
RUN_KEY_RENAME_MAP = {"sampleId": "sample_accession", "id": "run_accession", "instrumentModel": "instrument_model"}
SAMPLE_ATTRIBUTE_MAP = {
    "geographic location (country and/or sea)": "country",
    "geographic location (latitude)": "latitude",
    "geographic location (longitude)": "longitude",
    "collection date": "collection_date",
}

RETRY_COUNT = 3

ENA_SUBMISSION_URL = "https://www.ebi.ac.uk/ena/submit/report"
ENA_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_BROWSER_URL = "https://www.ebi.ac.uk/ena/browser/api/xml"

RUN_DEFAULT_FIELDS = ",".join(["secondary_study_accession", "sample_accession"])

STUDY_RUN_DEFAULT_FIELDS = ",".join(["sample_accession", "run_accession", "instrument_model"])

ASSEMBLY_DEFAULT_FIELDS = "sample_accession"

SAMPLE_DEFAULT_FIELDS = ",".join(["sample_accession", "collection_date", "country", "location"])

STUDY_DEFAULT_FIELDS = "study_accession,study_description"
