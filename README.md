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
IF you don't wish to involve conda at all, due to e.g. server restrictions or hosting on third party platforms (such as Galaxy), you could look through the `serovar_detector/workflow/envs` files and try to replicate an environment containing these tools as well as the tools listed in **Requirements** section (conda is only required when letting Snakemake handle the pipeline environments). In theory if these tool can exist within the same environment, you can skip the Snakemake conda handling **entirely** by adding the `-c` flag (e.g. `serovar_detector ... -c`), it's untested and unsupported though.

# Usage
Serovar detector runs with different types of modes.
* Single sample reads (`-1` and `-2`)
* Single sample assembly (`-A`)
* Batch samples reads (`-r`)
* Batch samples assemblies (`-a`)

## Input
Serovar detector assumes use of **either**/**and** Illumina Paired end **gz**-compressed sequencing data _(_R1|R2.fastq.gz)_ **or**/**and** preassembled genomes _(.fasta)_.

### Batch mode
Serovar detector will scan the top folder of the specified reads directory and assemblies directory when using the `-r` and `-a` options respectively.
```
serovar_detector -r /path/to/input/reads -a /path/to/input/assemblies -o /path/to/output -t 3
```

### Single sample mode
When running single sample mode (e.g. for setting up modules on galaxy or integrating Serovar detector in other pipelines) you can specify the file path for both read mates using `-1` and `-2` options, and/or the file path for the preassembled genome using the `-A` option.
```
serovar_detector -1 /path/to/input/sample_R1.fastq.gz -2 /path/to/input/sample_R2.fastq.gz -A /path/to/input/sample_assembly.fasta -o /path/to/output -t 3
```

### Mix of both
It's possible to mix batch and single mode, e.g. by specifying the `-r` (reads_dir) option and `-A` (assembly_file) option. Then all samples detected in the read directory alongside the single assembly will be included in a given run.
```
serovar_detector -r /path/to/input/reads -A /path/to/input/sample_assembly.fasta -o /path/to/output -t 3
```
If specifying both single sample AND batch type of the same data type (in this example assemblies) the batch mode will be ignored (.e.g specifying `-A` will cause `-a` to be ignored).

## Output
Serovar detector uses [KMerAligner (KMA)](https://bitbucket.org/genomicepidemiology/kma) to map raw reads against the capsule gene database of Serovar detector, and it uses [Blastn](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) to map the capsule gene database against the assembled genomes. Please note that it is possible to include both raww reads and assemblies of the same sample, this would lead to both mapping results being reported separately.

A summary of detected capsule genes and their derived serovars are provided for each individual _sample_ and _mapper_, in a single tab-separated file: `/path/to/output/serovars.tsv`.

## Options
```
serovar_detector -h
usage: serovar_detector [-h] [-1 --r1] [-2 --r2] [-A --assembly]
                        [-r --reads_dir] [-a --assembly_dir] -o --outdir
                        [-T --theshold] [-t --threads] [-k] [-F] [-c] [-n]
                        [-D --database] [-d]

Screen read files and assemblies for Serovar biomarker genes, in order to
preovide suggestions for isolate serovar. Currently only supporting
Actinobacillus Pleuropneumoniae.

options:
  -h, --help         show this help message and exit
  -1 --r1            Path to sample read mate 1 (Disables --reads_dir)
  -2 --r2            Path to sample read mate 2 (Disables --reads_dir)
  -A --assembly      Path to sample assembly (Disables --assembly_dir)
  -r --reads_dir     Input path to reads directory
  -a --assembly_dir  Input path to assembly directory
  -o --outdir        Output path to Results and Temporary files directory
  -T --theshold      Cutoff threshold of match coverage and identity. Ignore
                     threshold by setting to 0 or False. (Default: 98)
  -t --threads       Number of threads to allocate for the pipeline. (Default:
                     3)
  -k                 Preserve temporary files such as KMA result files.
                     (Default: False)
  -F                 Force rerun of all tasks in pipeline. (Default: False)
  -c                 Don't let snakemake handle conda execution in rules.
                     Enable this option if the pipeline should run in the
                     current loaded environment. (Default: False)
  -n                 Perform a dry run with Snakemake to see jobs but without
                     executing them. (Default: False)
  -D --database      Path and prefix to kmer-aligner database. (Default: /home
                     /cucumbergebt/micromamba/envs/serovar_detector/lib/python
                     3.14/site-packages/serovar_detector/db/Actinobacillus_ple
                     uropneumoniae)
  -d                 Enable debug mode, prints more messages and stores
                     snakemake object for inspection in R. (Default: False)
```

# Issues or questions
If you encounter any issues or have any questions, you are more than welcome to post these in the issues section of this repository (Requires a GitHub account).

# Citation
If you are using our tool in your analysis, please consider to cite us.

> Angen Ø, Karstensen KT, Vilaró A, et al. Serotyping of *Actinobacillus pleuropneumoniae* based on whole genome sequencing: validation of a bioinformatic tool. *Microb Genom.* 2025;11(7):001434. doi:10.1099/mgen.0.001434
