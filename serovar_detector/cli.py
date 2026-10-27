#!/bin/python

from importlib import resources
from pathlib import Path
import argparse
from .__version__ import __version__
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
  print("Generating config file", end = "...")
  # Define config file
  config_file = outdir / "config.yaml"

  # Check for existing config file
  if not config_file.is_file():
    config_dir = config_file.parent

    # Check for existing config dir
    if not config_dir.is_dir():
      print("No config dir detected, creating directory.")
      config_dir.mkdir(parents = True)

  out_path = outdir.absolute()

  config = {
    "database": str(database), "outdir": str(out_path), "threshold": threshold,
    "tmpdir": str(tmpdir), "debug": debug, "threads": threads
  }

  with open(config_file, "w") as config_yaml:
    yaml.dump(config, config_yaml)

  if config_file.exists():
    print(f"Success: Overwritting {config_file}")
  else:
    print(f"Success: Written to {config_file}")

  return config_file


def create_file_table(pattern, file_list, type):
  # Search file names for metadata
  search = [re.search(pattern, str(file_name)) for file_name in file_list]

  # Generate DataFrame from search object groups
  metadata_raw = pandas.DataFrame({
    "sample_name": [search.group("sample_name") for search in search],
    "file": [search.group("file") for search in search],
    "type": type
  })


  if type == "Reads":
    metadata_raw["mate"] = [search.group("mate") for search in search]
    metadata = metadata_raw.sort_values(by = ["sample_name", "mate"])
  else:
    metadata = metadata_raw.sort_values(by = "sample_name")

  # Sort table
  

  return metadata




def screen_files(path, type, single):
  read_pattern = r"^(?P<file>\S+\/(?P<sample_name>\S+?)((_S\d+)(_L\d+))?_(?P<mate>[Rr]?[12])(_\d{3})?\.(?P<ext>((fastq)?(fq)?)?(\.gz)?))"
  assembly_pattern = r"^(?P<file>\S+\/(?P<sample_name>\S+)\.(?P<ext>[fast]+))"

  if type == "Reads" and not single:
    # Screen for files
    fastqs = [sample_read for sample_read in glob.glob(f"{path}/*.fastq*")]
    fqs = [sample_read for sample_read in glob.glob(f"{path}/*.fq*")]

    # Combine file lists
    read_files = fastqs + fqs

    # Create metadata from read_dir
    metadata = create_file_table(pattern = read_pattern, file_list = read_files, type = type)

  elif type == "Reads" and path: 
    # Create metadata from single read pair
    metadata = create_file_table(pattern = read_pattern, file_list = path, type = type)

  elif type == "Assembly" and not single:
    # Screen for files
    fastas = [sample_assembly for sample_assembly in glob.glob(f"{path}/*.fasta")]
    fas = [sample_assembly for sample_assembly in glob.glob(f"{path}/*.fa")]

    # Combine file lists
    assembly_files = fastas + fas

    # Create metadata from assembly_dir
    metadata = create_file_table(pattern = assembly_pattern, file_list = assembly_files, type = type)

  elif type == "Assembly" and path:
    # Create metadata from single assembly
    metadata = create_file_table(pattern = assembly_pattern, file_list = path, type = type)

  else:
    # Create empty metadata
    metadata = pandas.DataFrame(columns = ["sample_name", "mate", "file", "type"])

  return metadata


def create_symlinks(metadata, tmpdir):
  print(f"Linking sample files", end = "... ")
  # Iterating over each row
  for row in range(len(metadata.index)):

    # Extracting sample information
    sample_metadata = metadata.loc[row]
    sample_name = sample_metadata["sample_name"]

    # Generate link filenames
    if sample_metadata["type"] == "Reads":
      read_dir = tmpdir / "reads"
      if not read_dir.is_dir():
        # Ensuring reads directory  exists
        read_dir.mkdir(parents = True)

      # Defining input and output
      mate = sample_metadata["mate"].replace("R", "")
      sample_file = Path(sample_metadata["file"]).resolve()
      sample_link = read_dir / f"{sample_name}_R{mate}.fastq.gz"

    elif sample_metadata["type"] == "Assembly":
      ass_dir = tmpdir / "assemblies"
      
      if not ass_dir.is_dir():
        # Ensuring assemblies directory exists
        ass_dir.mkdir(parents = True)

      # Defining input and output
      sample_file = Path(sample_metadata["file"]).resolve()
      sample_link = ass_dir / f"{sample_name}.fasta"

    # Symlinking file
    if not sample_link.exists():
      try:
        sample_link.symlink_to(sample_file)
      except FileExistsError:
        sample_link.unlink()
        sample_link.symlink_to(sample_file)

  print(f"Sucess: Files successfully linked to {tmpdir}!")
  return True


def make_sample_sheet(table):
  print("Generating sample sheet", end = "... ")

  file_count = len(table.index)

  if file_count > 0:
    sample_subset = table[["sample_name"]]
    sample_sheet = sample_subset.drop_duplicates()
  else:
    print("Failed: Subsample sheet has no rows")
    sample_sheet = False

  print(f"Success: A total of {len(sample_sheet.index)} samples have been annotated!")
  return sample_sheet


def write_sample_sheet(sample_sheet, pepdir):
  sample_file = pepdir / "sample_sheet.csv"

  print("Writting sample sheet", end = "... ")
  sample_sheet.to_csv(sample_file, index = False)

  if sample_file.exists():
    print(f"Success: Overwritting {sample_file}")
  else:
    print(f"Success: Written to {sample_file}")


def write_subsample_sheet(subsample_sheet, pepdir):
  subsample_file = pepdir / "subsample_sheet.csv"

  print("Writting subsample sheet", end = "... ")
  subsample_sheet.to_csv(subsample_file, index = False)

  if subsample_file.exists():
    print(f"Success: Overwritting {subsample_file}")
  else:
    print(f"Success: Written to {subsample_file}")


def write_PEP(pepdir):
  # Define pep files
  pep_file = pepdir / "project_config.yaml"

  # Generate PEP configuration:
  PEP_header = "pep_version: 2.1.0\n"
  PEP_sample = "sample_table: 'sample_sheet.csv'\n"
  PEP_subsample = "subsample_table: 'subsample_sheet.csv'"

  print("Writing project configuration file", end = "... ")
  with open(pep_file, "w") as config_file:
    config_file.write(PEP_header)
    config_file.write(PEP_sample)
    config_file.write(PEP_subsample)

  if pep_file.exists():
    print(f"Success: Overwriting {pep_file}")
  else:
    print(f"Success: Written to {pep_file}")


def generate_sheets(r1, r2, reads_dir, assembly, assembly_dir, outdir, tmpdir):
  
  if not outdir.is_dir():
    # Generating dirs
    outdir.mkdir(parents = True)
  if not tmpdir.is_dir():
    tmpdir.mkdir(parents = True)

  # Determining single vs batch
  reads_single_sample = r1 and r2

  if reads_single_sample:
    reads_metadata = screen_files(path = [r1, r2], type = "Reads", single = True)
  else:
    reads_metadata = screen_files(path = reads_dir, type = "Reads", single = False)
  
  if assembly:
    assembly = [assembly]
    assembly_metadata = screen_files(path = assembly, type = "Assembly", single = True)
  else:
    assembly_metadata = screen_files(path = assembly_dir, type = "Assembly", single = False)


  metadata = pandas.concat([reads_metadata, assembly_metadata], ignore_index = True, sort = True)

  # Ensuring not all samples have been filtered out
  sample_size = len(metadata.index)
  if sample_size > 0:

    # Generate symlinks
    create_symlinks(metadata, tmpdir)

    # Ensure PEP folder exists
    pepdir = Path(outdir / "schemas")

    if not pepdir.is_dir():
      pepdir.mkdir()


    sample_sheet = make_sample_sheet(metadata)
    sample_sheet_updated = write_sample_sheet(sample_sheet, pepdir)

    subsample_sheet = metadata[["sample_name", "type", "mate", "file"]]
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
  parser.add_argument("-r", "--reads_dir", metavar = "--reads_dir", dest = "reads_dir", help = "Input path to reads directory", required = False)
  parser.add_argument("-1", "--r1", metavar = "--r1", dest = "r1", help = "Path to sample read mate 1 (Disables --reads_dir)", required = False, default = None)
  parser.add_argument("-2", "--r2", metavar = "--r2", dest = "r2", help = "Path to sample read mate 2 (Disables --reads_dir)", required = False, default = None)
  parser.add_argument("-a", "--assembly_dir", metavar = "--assembly_dir", dest = "assembly_dir", help = "Input path to assembly directory", required = False)
  parser.add_argument("-A", "--assembly", metavar = "--assembly", dest = "assembly", help = "Path to sample assembly (Disables --assembly_dir)", required = False, default = None)
  parser.add_argument("-o", "--outdir", metavar = "--outdir", dest = "outdir", help = "Output path to Results and Temporary files directory", required = True)
  parser.add_argument("-T", "--theshold", metavar = "--theshold", dest = "threshold", help = "Cutoff threshold of match coverage and identity. Ignore threshold by setting to 0 or False. (Default: %(default)s)", default = 98)
  parser.add_argument("-t", "--threads", metavar = "--threads", dest = "threads", help = "Number of threads to allocate for the pipeline. (Default: %(default)s)", default = 3)
  parser.add_argument("-k", "--keep_tmp", dest = "keep_tmp", help = "Preserve temporary files such as KMA result files. (Default: %(default)s)", action = "store_true")
  parser.add_argument("-F", "--force", dest = "force", help = "Force rerun of all tasks in pipeline. (Default: %(default)s)", action = "store_true")
  parser.add_argument("-c", "--noconda", dest = "noconda", help = "Don't let snakemake handle conda execution in rules. Enable this option if the pipeline should run in the current loaded environment. (Default: %(default)s)", action = "store_true")
  parser.add_argument("-n", "--dry_run", dest = "dry_run", help = "Perform a dry run with Snakemake to see jobs but without executing them. (Default: %(default)s)", action = "store_true")
  parser.add_argument("-D", "--database", metavar = "--database", dest = "database", help = "Path to database Fasta file. (Default: %(default)s)", default = f"{PKG_DIR}/db/Actinobacillus_pleuropneumoniae.fasta")
  parser.add_argument("-d", "--debug", dest = "debug", help = "Enable debug mode, prints more messages and stores snakemake object for inspection in R. (Default: %(default)s)", action = "store_true")
  parser.add_argument("-v", "--version", action="version", version=f"serovar_detector {__version__}")

  return(parser.parse_args())


def main():

  # Derrive arguments
  args = parse_arguments()

  database = Path(args.database).absolute()
  outdir = Path(args.outdir).absolute()
  threshold = args.threshold
  keep_tmp = args.keep_tmp
  threads = int(args.threads)
  force = args.force
  noconda = args.noconda
  dry_run = args.dry_run
  debug = args.debug
  tmpdir = Path(outdir / "tmp")


  # Polish input
  r1 = None
  r2 = None
  reads_dir = None
  if args.r1:
    r1 = Path(args.r1).absolute()
    r2 = Path(args.r2).absolute()
  elif args.reads_dir:
      reads_dir = Path(args.reads_dir).absolute()

  assembly = None
  assembly_dir = None

  if args.assembly:
    assembly = Path(args.assembly).absolute()
  elif args.assembly_dir:
      assembly_dir = Path(args.assembly_dir).absolute()

  if not (r1 or r2 or reads_dir or assembly or assembly_dir):
    raise ValueError("No input files or directories provided for reads or assemblies.")

  # Prepare config file for snakemake
  config_file = generate_configfile(database = database, outdir = outdir, threshold = threshold, threads = threads, debug = debug, tmpdir = tmpdir)

  # Generate subsample sheet
  sample_files = generate_sheets(r1 = r1, r2 = r2, reads_dir = reads_dir, assembly = assembly, assembly_dir = assembly_dir, outdir = outdir, tmpdir = tmpdir)


  if len(sample_files) == 0:
    print("Nothing to do exitting!")
    sys.exit(0)

  snake_args = f"--snakefile {PKG_DIR}/workflow/Snakefile --configfile {config_file} --cores {threads} "
  if not noconda:
    snake_args += f"--use-conda --conda-prefix {PKG_DIR}/conda "
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
    sys.exit(1)
