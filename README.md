# Serovar detector
Actinobacillus pleuropneumoniae causes severe respiratory illness in production pigs and piglets. One major task in the prohobition of spread as well as treatment, is to identify the composition of capsule genes. These capsule genes can in combination be used too provide servariant typing of A. pleuropneumoniae.

This repository provides a pipeline for detecting capsule genes and dessiminate serovar from combination of present genes.

# Setup
## Requirements
* snakemake
* conda >= 24.7.1 

## Installation
1. (Recommendation) Use Micromamba or the like to set up an environment with the above requirement.
2. Clone repository to desired location (e.g. ~/repos)
```
git clone https://github.com/KasperThystrup/serovar_detector.git ~/repos/serovar_detector
```

# Usage
## Quick start
Assuming that you have navigated into the repository folder (e.g. ~/repos/serovar_detector) and are executing the command from an conda environment with installed requirements.
```
python serovar_detector.py -r /path/to/input/reads -a /path/to/input/assemblies -D db/Actinobacillus_pleuropneumoniae -o /path/to/output/serovar_detector/ -t [nr. of threads]
```

## Options
```
usage: serovar_detector.py [-h] [-r --reads_dir] [-a --assembly_dir] -D
                           --database -o --outdir [-T --theshold] [-R] [-b]
                           [-B] [-k] [-t --threads] [-F] [-n] [-d]

Screen read files and assemblies for Serovar biomarker genes, in order to
preovide suggestions for isolate serovar. Currently only supporting
Actinobacillus Pleuropneumoniae.
```

An overview of the available options:
```
  -h, --help         show this help message and exit
  -r --reads_dir     Input path to reads directory
  -a --assembly_dir  Input path to assembly directory
  -D --database      Path and prefix to kmer-aligner database
  -o --outdir        Output path to Results and Temporary files directory
  -T --theshold      Cutoff threshold of match coverage and identity. Ignore
                     threshold by setting to 0 or False. (Default 98)
  -R                 Append to existing results file. (Default False)
  -b                 Update existing blacklist file with new samples. Creates
                     a blacklist file if non exists. (Default False)
  -B                 Ignore and overwrite existing blacklist file. Creates a
                     blacklist if non exists. (Default False)
  -k                 Preserve temporary files such as KMA result files.
                     (Default False)
  -t --threads       Number of threads to allocate for the pipeline. (Default
                     3)
  -F                 Force rerun of all tasks in pipeline. (Default False)
  -n                 Perform a dry run with Snakemake to see jobs but without
                     executing them. (Default False)
  -d                 Enable debug mode, stores snakemake object for inspection
                     in R. (Default False)
```

