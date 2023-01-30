# Public bins and MAGs uploader
Python script to prepare bins and MAGs upload to ENA in fasta format (European Nucleotide Archive). This script generates xmls and manifests necessary for submission with webin-cli. 

It takes as input one tsv (tab-separated values) table expecting the following columns:
  * _genome_name_: genome id (unique string identifier, shorter than 20 characters)
  * _run_accessions_: run(s) genome was generated from (DRR/ERR/SRRxxxxxx accessions)
  * _assembly_software_: assemblerName_vX.X
  * _binning_software_: binnerName_vX.X
  * _binning_parameters_: binning parameters
  * _stats_generation_software_: software_vX.X
  * _completeness_: `float`
  * _contamination_: `float`
  * _rRNA_presence_: `True/False` if 5S, 16S, and 23S genes have been detected in the genome
  * _NCBI_lineage_: full NCBI lineage, either in tax ids (`integers`) or `strings`. Format: x;y;z;...
  * _metagenome_: needs to be listed in the taxonomy tree [here](<https://www.ebi.ac.uk/ena/browser/view/408169?show=tax-tree>) (you might need to press "Tax tree - Show" in the right most section of the page)
  * _co_-assembly: `True/False`, whether the genome was generated from a co-assembly
  * _genome_coverage_ : genome coverage
  * _genome_path_: path to genome to upload (already compressed)
  * _broad_environment_: `string` (explanation following)
  * _local_environment_: `string` (explanation following)
  * _environmental_medium_: `string` (explanation following)

According to ENA checklist's guidelines, 'broad_environment' describes the broad ecological context of a sample - desert, taiga, coral reef, ... 'local_environment' is more local - lake, harbour, cliff, ... 'environmental_medium' is either the material displaced by the sample, or the one in which the sample was embedded prior to the sampling event - air, soil, water, ... 
For host-associated metagenomic samples, variables can be defined similarly to the following example for the chicken gut metagenome: 'Biome: chicken digestive system, Feature: digestive tube, Material: caecum. More information can be found at <https://www.ebi.ac.uk/ena/browser/view/ERC000050> for bins and [ERC000047](<https://www.ebi.ac.uk/ena/browser/view/ERC000047>) for MAGs under field names "broad-scale environmental context", "local environmental context", "environmental medium"

For the script to work, raw-read runs should already be available on the INSDC (ENA, NCBI, or DDBJ), hence at least one DRR|ERR|SRR accession should be available for every genome to be uploaded. 

Files to be uploaded will need to be compressed (e.g. already in .gz format). 

## How to run
The script needs `python`, `pandas`, and `requests` to run. A quick way of creating an environment is via [`venv`](<https://virtualenv.pypa.io/en/latest/installation//>) (e.g. via `apt install python3-virtualenv`):
```bash
# Create python environment
virtualenv -p python3 venv

# Source environment and install requirements
source venv/bin/activate && pip install -r requirements.txt
```

After this, you just need to run the script as follows:

```bash
python genome_upload.py -u UPLOAD_STUDY --genome_info METADATA_FILE (--mags | --bins) --xmls --manifests --webin WEBIN_ID --password PASSWORD [--out] [--force] [--live]
```

where

  * `-u UPLOAD_STUDY`: study accession for genomes upload to ENA (in format ERPxxxxxx or PRJEBxxxxxx)
  * `---genome_info METADATA_FILE` : genomes metadata file in tsv format
  * `-m, --mags, --b, --bins`: select for bin or MAG upload
  * `--xmls`: creates submission and genome registration xmls
  * `--manifests`: creates a manifest file for every genome to upload (if xmls have already been generated and do not need to be updated, `--manifests` can be used without `--xmls`)
  * `--out`: output folder (default: working directory)
  * `--force`: forces reset of sample xml's backups
  * `--live`: registers genomes on ENA's live server. Omitting this option allows to validate samples beforehand
  * `--webin WEBIN_ID`: webin id (format: Webin_XXXXX)
  * `--password PASSWORD`: webin password
  * `--centre_name CENTRE_NAME`: name of the centre uploading genomes

### Produced files:
The script produces the following files and folders:
  * `bin_upload` or `MAG_upload` folder (according to the upload type) containing:
     - `manifests` folder: contains all generated manifests
     - `genome_samples.xml`: xml generated to register samples on ENA before the upload
     - `ENA_backup.json`: backup file to prevent re-download of metadata from ENA. Deletion can be forced with the `--force` option
     - `registered_bins\MAGs.tsv`: contains a list of genomes registered on ENA. This file is needed for manifest generation - do not delete it. If the submission hasn't been launched in `--live` mode, a `test` file with test accessions will be generated. 
     - `submission.xml`: xml needed for genome registration on ENA

## What to do next
Once manifest files are generated, it will be necessary to use ENA's webin-cli resource to upload genomes. More information can be found [here](<https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html>). 