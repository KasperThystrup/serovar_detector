#!/bin/python

from importlib import resources
import argparse
import os
import sys
import yaml
import glob
import subprocess
import re
import pandas
import shutil

# Define package installation folder

PKG_DIR = resources.files("serovar_detector")

def generate_configfile(database, outdir, threshold, threads, debug, tmpdir):
  # Define config file
  config_file = f"{outdir}/config.yaml"

  # Check for existing config file
  if not os.path.isfile(config_file):
    config_dir = os.path.dirname(config_file)

    # Check for existing config dir
    if not os.path.isdir(config_dir):
      print("No config dir detected, creating directory.")
      os.makedirs(config_dir)

  out_path = os.path.abspath(outdir).rstrip("/")

  config = {
    "database": database, "outdir": out_path, "threshold": threshold,
    "tmpdir": tmpdir, "debug": debug, "threads": threads
  }

  with open(config_file, "w") as config_yaml:
    yaml.dump(config, config_yaml)

  return config_file


def screen_files(directory, type):
  if type == "Reads":
    # Screen for files
    fastqs = [sample_read for sample_read in glob.glob(f"{directory}/*.fastq*")]
    fqs = [sample_read for sample_read in glob.glob(f"{directory}/*.fq*")]

    # Combine file lists
    read_files = fastqs + fqs

    # Search file names for metadata
    pattern = r"^(?P<file>\S+\/(?P<sample_name>\S+?)((_S\d+)(_L\d+))?_(?P<mate>[Rr]?[12])(_\d{3})?\.(?P<ext>((fastq)?(fq)?)?(\.gz)?))"
    search = [re.search(pattern, read_file) for read_file in read_files]

    # Generate DataFrame from search object groups
    metadata_raw = pandas.DataFrame({
      "sample_name": [search.group("sample_name") for search in search],
      "mate": [search.group("mate") for search in search],
      "file": [search.group("file") for search in search],
      "type": type
    })

    # Sort table
    metadata = metadata_raw.sort_values(by = ["sample_name", "mate"])

  elif type == "Assembly":
    # Screen for files
    fastas = [sample_assembly for sample_assembly in glob.glob(f"{directory}/*.fasta")]
    fas = [sample_assembly for sample_assembly in glob.glob(f"{directory}/*.fa")]

    # Combine file lists
    assembly_files = fastas + fas

    # Search file names for metadata
    pattern = r"^(?P<file>\S+\/(?P<sample_name>\S+)\.(?P<ext>[fast]+))"
    search = [re.search(pattern, assembly_file) for assembly_file in assembly_files]

    # Generate DataFrame from search object groups
    metadata_raw = pandas.DataFrame({
      "sample_name": [search.group("sample_name") for search in search],
      "file": [search.group("file") for search in search],
      "type": type
    })

    # Sort table
    metadata = metadata_raw.sort_values(by = "sample_name")

  else:
    return pandas.DataFrame({"sample_name", "mate", "file", "type"})

  return metadata #file_metadata


def create_symlinks(metadata, outdir):
  # Iterating voer each row
  for row in range(len(metadata.index)):
    # Extracting sample information
    sample_metadata = metadata.loc[row]
    sample_name = sample_metadata["sample_name"]

    # Generate link filenames
    if sample_metadata["type"] == "Reads":
      # Ensuring outdir exists
      os.makedirs(name = f"{outdir}/reads", exist_ok = True)

      # Defining input and output
      mate = sample_metadata["mate"]
      sample_file = os.path.realpath(sample_metadata["file"])
      sample_link = f"{outdir}/reads/{sample_name}_{mate}.fastq.gz"

    else:
      # Ensuring outdir exists
      os.makedirs(name = f"{outdir}/assemblies", exist_ok = True)

      # Defining input and output
      sample_file = os.path.realpath(sample_metadata["file"])
      sample_link = f"{outdir}/assemblies/{sample_name}.fasta"

    # Symlinking file
    if not os.path.exists(sample_link):
      try:
        os.symlink(src = sample_file, dst = sample_link)
      except FileExistsError:
        os.unlink(sample_link)
        os.symlink(src = sample_file, dst = sample_link)

  print("Files successfully linked!")
  return True


def make_sample_sheet(table):
  print("Generating sample sheet", end = "... ")

  file_count = len(table.index)

  if file_count > 0:
    sample_subset = table[["sample_name", "type"]]
    sample_sheet = sample_subset.drop_duplicates()
  else:
    print("Failed: Subsample sheet has no rows")
    sample_sheet = False

  print(f"Success: A total of {len(sample_sheet.index)} samples have been annotated!")
  return sample_sheet


def write_sample_sheet(sample_sheet, pepdir):
  sample_file = f"{pepdir}/sample_sheet.csv"
  sample_exists = os.path.exists(sample_file)

  print("Writting sample sheet", end = "... ")
  sample_sheet.to_csv(sample_file, index = False)
  if sample_exists:
    print(f"Success: Overwritting {sample_file}")
  else:
    print(f"Success: Written to {sample_file}")


def write_subsample_sheet(subsample_sheet, pepdir):
  subsample_file = f"{pepdir}/subsample_sheet.csv"
  subsample_exists = os.path.exists(subsample_file)

  print("Writting subsample sheet", end = "... ")
  subsample_sheet.to_csv(subsample_file, index = False)

  if subsample_exists:
    print(f"Success: Overwritting {subsample_file}")
  else:
    print(f"Success: Written to {subsample_file}")


def write_PEP(pepdir):
  # Define pep files
  pep_file = f"{pepdir}/project_config.yaml"
  pep_exists = os.path.exists(pep_file)

  # Generate PEP configuration:
  PEP_header = "pep_version: 2.1.0\n"
  PEP_sample = "sample_table: 'sample_sheet.csv'\n"
  PEP_subsample = "subsample_table: 'subsample_sheet.csv'"

  print("Writing project configuration file", end = "... ")
  with open(pep_file, "w") as config_file:
    config_file.write(PEP_header)
    config_file.write(PEP_sample)
    config_file.write(PEP_subsample)

  if not pep_exists:
    print(f"Success: Written to {pep_file}")
  elif pep_exists:
    print(f"Success: Overwriting {pep_file}")


def generate_sheets(reads_dir, assembly_dir, outdir, tmpdir):
  # Generating dirs
  os.makedirs(outdir, exist_ok=True)
  os.makedirs(tmpdir, exist_ok=True)

  # Screening input files
  reads_metadata = screen_files(directory = reads_dir, type = "Reads")
  assembly_metadata = screen_files(directory = assembly_dir, type = "Assembly")

  metadata = pandas.concat([reads_metadata, assembly_metadata], ignore_index = True, sort = True)

  # Ensuring not all samples have been filtered out
  sample_size = len(metadata.index)
  if sample_size > 0:

    # Generate symlinks
    create_symlinks(metadata, tmpdir)

    # Ensure PEP folder exists
    pepdir = f"{outdir}/schemas"
    pepdir_exists = os.path.isdir(pepdir)

    if not pepdir_exists:
      os.makedirs(pepdir, exist_ok = True)

    sample_sheet = make_sample_sheet(metadata)
    sample_sheet_updated = write_sample_sheet(sample_sheet, pepdir)

    subsample_sheet = metadata[["sample_name", "mate", "file"]]
    subsample_sheet_updated = write_subsample_sheet(subsample_sheet, pepdir)
    sample_files = metadata['file']

    # Generate PEP configuration files.
    write_PEP(pepdir)
  else:
    print("No samples detected.")
    sample_files = []

  return(sample_files)


def parse_arguments():
  parser = argparse.ArgumentParser(description = "Screen read files and assemblies for Serovar biomarker genes, in order to preovide suggestions for isolate serovar. Currently only supporting Actinobacillus Pleuropneumoniae.")
  parser.add_argument("-r", metavar = "--reads_dir", dest = "reads_dir", help = "Input path to reads directory", required = False)
  parser.add_argument("-a", metavar = "--assembly_dir", dest = "assembly_dir", help = "Input path to assembly directory", required = False)
  parser.add_argument("-D", metavar = "--database", dest = "database", help = "Path and prefix to kmer-aligner database. (Default: %(default)s)", default = f"{PKG_DIR}/db/Actinobacillus_pleuropneumoniae")
  parser.add_argument("-o", metavar = "--outdir", dest = "outdir", help = "Output path to Results and Temporary files directory", required = True)
  parser.add_argument("-T", metavar = "--theshold", dest = "threshold", help = "Cutoff threshold of match coverage and identity. Ignore threshold by setting to 0 or False. (Default: %(default)s)", default = 98)
  parser.add_argument("-t", metavar = "--threads", dest = "threads", help = "Number of threads to allocate for the pipeline. (Default: %(default)s)", default = 3)
  parser.add_argument("-k", dest = "keep_tmp", help = "Preserve temporary files such as KMA result files. (Default: %(default)s)", action = "store_true")
  parser.add_argument("-F", dest = "force", help = "Force rerun of all tasks in pipeline. (Default: %(default)s)", action = "store_true")
  parser.add_argument("-c", dest = "noconda", help = "Don't let snakemake handle conda execution in rules. Enable this option if the pipeline should run in the current loaded environment. (Default: %(default)s)", action = "store_true")
  parser.add_argument("-n", dest = "dry_run", help = "Perform a dry run with Snakemake to see jobs but without executing them. (Default: %(default)s)", action = "store_true")
  parser.add_argument("-d", dest = "debug", help = "Enable debug mode, prints more messages and stores snakemake object for inspection in R. (Default: %(default)s)", action = "store_true")

  return(parser.parse_args())


def main():

  # Derrive arguments
  args = parse_arguments()

  database = os.path.abspath(args.database)
  outdir = os.path.abspath(args.outdir)
  threshold = args.threshold
  keep_tmp = args.keep_tmp
  threads = int(args.threads)
  force = args.force
  noconda = args.noconda
  dry_run = args.dry_run
  debug = args.debug
  tmpdir = f"{outdir}/tmp"


  # Polish input
  if not args.reads_dir:
    reads_dir = args.reads_dir
  else:
    reads_dir = os.path.abspath(args.reads_dir)
  if not args.assembly_dir:
    assembly_dir = args.assembly_dir
  else:
    assembly_dir = os.path.abspath(args.assembly_dir)

  # Prepare config file for snakemake
  config_file = generate_configfile(database = database, outdir = outdir, threshold = threshold, threads = threads, debug = debug, tmpdir = tmpdir)

  # Generate subsample sheet
  sample_files = generate_sheets(reads_dir = reads_dir, assembly_dir = assembly_dir, outdir = outdir, tmpdir = tmpdir)


  if len(sample_files) == 0:
    print("Nothing to do exitting!")
    sys.exit(0)

  snake_args = f"--snakefile {PKG_DIR}/workflow/Snakefile --configfile {config_file} --cores {threads} "
  if not noconda:
    snake_args += "--use-conda "
  if force:
    snake_args += "--forceall "
  if dry_run:
    snake_args += "--dryrun "


  snakemake_cmd = f"snakemake {snake_args}"
  if debug:
    print(f"Executing: {snakemake_cmd}")

  snake_success = subprocess.Popen(snakemake_cmd, shell = True).wait()

  if snake_success == 0:
    if not keep_tmp:
      print("Cleaning up temporary files.")
      shutil.rmtree(tmpdir)
    else:
      print(f"Keeping all temporary files, inspect them at: {tmpdir}")

    print("All Done!")
  else:
    print("Something went wrong while executing snakemake.")
