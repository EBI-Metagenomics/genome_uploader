# Public bins and MAGs uploader
Python script to prepare bins and MAGs upload to ENA (European Nucleotide Archive). This script generates xmls and manifests necessary for submission with webin-cli. 

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
  * _metagenome_: needs to be listed in the taxonomy tree at <https://www.ebi.ac.uk/ena/browser/view/408169?show=tax-tree>
  * _co_-assembly: `True/False`, whether the genome was generated from a co-assembly
  * _genome_coverage_ : genome coverage
  * _genome_path_: path to genome to upload
  * _broad_environment_: `string` (explanation following)
  * _local_environment_: `string` (explanation following)
  * _environmental_medium_: `string` (explanation following)

According to ENA checklist's guidelines, 'broad_environment' describes the broad ecological context of a sample - desert, taiga, coral reef, ... 'local_environment' is more local - lake, harbour, cliff, ... 'environmental_medium' is either the material displaced by the sample, or the one in which the sample was embedded prior to the sampling event - air, soil, water, ... For host-associated metagenomic samples, variables can be defined based on our biome tree. For example, for chicken gut metagenome: 'Biome: chicken digestive system, Feature: digestive tube, Material: caecum. More information can be found at <https://www.ebi.ac.uk/ena/browser/view/ERC000050> for bins and [ERC000047](<https://www.ebi.ac.uk/ena/browser/view/ERC000047) for MAGs under field names "broad-scale environmental context", "local environmental context", "environmental medium"

## How to execute
The script needs `python` to run.

```bash
python genome_upload.py -u UPLOAD_STUDY --genome_info METADATA_FILE (--mags | --bins) --xmls --manifests --webin --password [--out] [--force] [--live]
```

where

  * `-u UPLOAD_STUDY`: study accession for genomes upload
  * `---genome_info METADATA_FILE` : genomes metadata  in tsv format
  * `-m, --mags, --b, --bins`: select for bin or MAG upload
  * `--xmls`: creates submission and genome registration xmls
  * `--manifests`: creates a manifest file for every genome to upload
  * `--out`: output folder (default: working directory)
  * `--force`: forces reset of sample xml's backups
  * `--live`: registers genomes on ENA's live server. Omitting this option allows to validate samples beforehand
  * `--webin WEBIN`: webin id
  * `--password PASSWORD`: webin password
  * `--centre_name CENTRE_NAME`: name of the centre uploading genomes