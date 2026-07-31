# Serovar detector
Actinobacillus pleuropneumoniae causes severe respiratory illness in production pigs and piglets. One major task in the prohobition of spread as well as treatment, is to identify the composition of capsule genes. These capsule genes can in combination be used too provide servariant typing of A. pleuropneumoniae.

This repository provides a pipeline for detecting capsule genes and dessiminate serovar from combination of present genes. 

# Quick start
```
# Install once only
micromamba create -n serovar_detector bioconda::serovar_detector --yes

# Execution
micromamba run -n serovar_detector serovar_detector -r /path/to/input/reads -a /path/to/input/assemblies -o /path/to/output
```
# Setup
## Requirements
While there are several other tools being used, standard usage assumes that Snakemake automatically should handle installation of downstream software.
Here's the requirement for Serovar detector:
* snakemake >= 8.+
* conda >= 24.7.1
* pandas
* peppy

## Micromamba/Conda
Serovar detector is hosted on bioconda, so requirements can be automatically fixed by using micromamba/conda. 
```
micromamba create -n serovar_detector bioconda::serovar_detector
```

## Pip installation
Installation through pip _can_ be a bit more tricky, it assumes you have conda preinstalled on your system (or in a virtual environment). Since Serovar detector is hosted on PyPi, it can be installed by simply running:
```
pip install serovar_detector
```

### Alternative to Conda
IF you don't wish to involve conda at all, due to e.g. server restrictions or hosting on third party platforms (such as Galaxy), you could look through the `serovar_detector/workflow/envs` files and try to replicate an environment containing these tools as well as the tools listed in **Requirements** section (except `conda`). In theory if these tool can exist within the same environment, you can skip conda handling entirely by adding the `-c` flag (e.g. `serovar_detector ... -c`), it's untested and unsupported though.

# Usage
```
serovar_detector -r /path/to/input/reads -a /path/to/input/assemblies -o /path/to/output -t 4
```
## Input
Serovar detector assumes use of either Illumina Paired end sequencing data or preassembled genomes. The top level of the `reads_dir/` and `assemblies_dir/` folders are screened for `.fastq.gz` and `.fasta` files respectively and sample name are derived automatically by ignoring the file extensions (and read mates).

## Output
Serovar detector uses [KMerAligner (KMA)](https://bitbucket.org/genomicepidemiology/kma) to map raw reads against the capsule gene database of Serovar detector, and it uses [Blastn](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) to map the capsule gene database against the assembled genomes. Please note that it is possible to include both raww reads and assemblies of the same sample, this would lead to both mapping results being reported separately.

A summary of detected capsule genes and their derived serovars are provided for each individual _sample_ and _mapper_, in a single tab-separated file: `/path/to/output/serovars.tsv`.

## Options
```
serovar_detector -h
usage: serovar_detector [-h] [-r --reads_dir] [-a --assembly_dir] [-D --database] -o --outdir [-T --theshold] [-t --threads] [-k] [-F] [-c] [-n] [-d]

Screen read files and assemblies for Serovar biomarker genes, in order to preovide suggestions for isolate serovar. Currently only supporting Actinobacillus Pleuropneumoniae.

options:
  -h, --help         show this help message and exit
  -r --reads_dir     Input path to reads directory
  -a --assembly_dir  Input path to assembly directory
  -D --database      Path and prefix to kmer-aligner database. (Default: /home/cucumbergebt/micromamba/envs/serovar_detector/lib/python3.14/site-
                     packages/serovar_detector/db/Actinobacillus_pleuropneumoniae)
  -o --outdir        Output path to Results and Temporary files directory
  -T --theshold      Cutoff threshold of match coverage and identity. Ignore threshold by setting to 0 or False. (Default: 98)
  -t --threads       Number of threads to allocate for the pipeline. (Default: 3)
  -k                 Preserve temporary files such as KMA result files. (Default: False)
  -F                 Force rerun of all tasks in pipeline. (Default: False)
  -c                 Don't let snakemake handle conda execution in rules. Enable this option if the pipeline should run in the current loaded environment. (Default: False)
  -n                 Perform a dry run with Snakemake to see jobs but without executing them. (Default: False)
  -d                 Enable debug mode, prints more messages and stores snakemake object for inspection in R. (Default: False)
```

# Issues or questions
If you encounter any issues or have any questions, you are more than welcome to post these in the issues section of this repository (Requires a GitHub account).

# Citation
If you are using our tool in your analysis, please consider to cite us.

> Angen Ø, Karstensen KT, Vilaró A, et al. Serotyping of *Actinobacillus pleuropneumoniae* based on whole genome sequencing: validation of a bioinformatic tool. *Microb Genom.* 2025;11(7):001434. doi:10.1099/mgen.0.001434
