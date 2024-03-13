import argparse
import os
import sys
import yaml
import glob
import subprocess
import re
import pandas
import shutil

def parse_arguments():
  parser = argparse.ArgumentParser(description = "Screen read files and assemblies for Serovar biomarker genes, in order to preovide suggestions for isolate serovar. Currently only supporting Actinobacillus Pleuropneumoniae.")
  parser.add_argument("-r", metavar = "--reads_dir", dest = "reads_dir", help = "Input path to reads directory", required = False)
  parser.add_argument("-a", metavar = "--assembly_dir", dest = "assembly_dir", help = "Input path to assembly directory", required = False)
  parser.add_argument("-D", metavar = "--database", dest = "database", help = "Path and prefix to kmer-aligner database", required = True)
  parser.add_argument("-T", metavar = "--theshold", dest = "threshold", help = "Cutoff threshold of match coverage and identity. Ignore threshold by setting to 0 or False. (Default 98)", default = 98)
  parser.add_argument("-o", metavar = "--outdir", dest = "outdir", help = "Output path to Results and Temporary files directory", required = True)
  parser.add_argument("-b", dest = "blacklisting", help = "Blacklist succesfsfully analysed samples, usable for surveillance / continous projects. (Default False)", action = "store_true")
  parser.add_argument("-B", dest = "clean_blacklisting", help = "Remove existing blacklist file. (Default False)", action = "store_true")
  parser.add_argument("-m", dest = "multithreading", help = "Enable multithreading during kmer alignment, use for huge samples only. (Default False)", action = "store_true")
  parser.add_argument("-k", dest = "keep_tmp", help = "Preserve temporary files such as KMA result files. (Default False)", action = "store_true")
  parser.add_argument("-c", dest = "keep_configs", help = "Preserve configuration files and sample_sheets. (Default False)", action = "store_true")
  parser.add_argument("-t", metavar = "--threads", dest = "threads", help = "Number of threads to allocate for the pipeline. (Default 3)", default = 3)
  parser.add_argument("-F", dest = "force", help = "Force rerun of all tasks in pipeline. (Default False)", action = "store_true")
  parser.add_argument("-f", dest = "force_results", help = "Force rerun of the Results task only. (Default False)", action = "store_true")
  parser.add_argument("-S", dest = "skipmake", help = "Skip Snakemake for requirering manual run of Snakemake. Config file will be generated.", action = "store_true")
  parser.add_argument("-n", dest = "dry_run", help = "Perform a dry run with Snakemake to see jobs but without executing them. (Default False)", action = "store_true")
  parser.add_argument("-d", dest = "debug", help = "Enable debug mode, stores snakemake object for inspection in R. (Default False)", action = "store_true")

  return(parser.parse_args())


def validate_snakemake(debug):
  here = os.listdir('.')
  workflow_here = 'workflow' in os.listdir('.')

  if (workflow_here):
    snakefile_here = "Snakefile" in os.listdir('workflow')

    if snakefile_here and debug:
      print("Snakefile detected")
    elif not snakefile_here:
      print('Error no Snakefile detected in workflow directory. Software is properbly corrupt, consider redownloading.')
      sys.exit(1)
  else:
    print('No workflow directory detected, are you sure you are running the script from the software folder?')
    sys.exit(1)


def generate_configfile(database, threshold, outdir, tmpdir, multithreading, threads, debug, skipmake):
  # Define config file
  config_file = "config/config.yaml"
  
  # Check for existing config file
  if not os.path.isfile(config_file):
    config_dir = os.path.dirname(config_file)
    
    # Check for existing config dir
    if not os.path.isdir(config_dir):
      print("No config dir detected, creating directory.")
      os.mkdir(config_dir)

  out_path = os.path.abspath(outdir).rstrip("/")
  
  config = {
    "database": database, "threshold": threshold,  "outdir": out_path, "tmpdir": tmpdir,
    "multithreading": multithreading, "debug": debug
  }

  with open(config_file, "w") as config_yaml:
    yaml.dump(config, config_yaml)


def screen_files(directory, type):
  if type == "Reads":
    # Screen for files
    fastqs = [sample_read for sample_read in glob.glob("%s/*.fastq*" %directory)]
    fqs = [sample_read for sample_read in glob.glob("%s/*.fq*" %directory)]

    # Combine file lists
    read_files = fastqs + fqs

    # Search file names for metadata
    pattern = "^(?P<file>\S+\/(?P<sample_name>\S+?)((_S\d+)?(_L\d+)?)?_(?P<mate>[Rr]?[12])(_\d{3})?\.(?P<ext>((fastq)?(fq)?)?(\.gz)?))"
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
    fastas = [sample_assembly for sample_assembly in glob.glob("%s/*.fasta" %directory)]
    fas = [sample_assembly for sample_assembly in glob.glob("%s/*.fa" %directory)]

    # Combine file lists
    assembly_files = fastas + fas

    # Search file names for metadata
    pattern = "^(?P<file>\S+\/(?P<sample_name>\S+)\.(?P<ext>[fast]+))"
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
    return(pandas.DataFrame({"sample_name", "mate", "file", "type"}))
  
  return(metadata) #file_metadata


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
      sample_file = sample_metadata["file"]
      sample_link = f"{outdir}/reads/{sample_name}_{mate}.fastq.gz"

    else:
      # Ensuring outdir exists
      os.makedirs(name = f"{outdir}/assemblies", exist_ok = True)

      # Defining input and output  
      sample_file = sample_metadata["file"]
      sample_link = f"{outdir}/assemblies/{sample_name}.fasta"

    # Symlinking file
    if not os.path.isfile(sample_link):
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

  print("Success: A total of %s samples have been annotated!" %len(sample_sheet.index))
  return(sample_sheet)


def write_sample_sheet(sample_sheet, pepdir, keep):
  sample_file = f"{pepdir}/sample_sheet.csv"
  sample_exists = os.path.exists(sample_file)

  print("Writting sample sheet", end = "... ")
  if not sample_exists or not keep:
    sample_sheet.to_csv(sample_file, index = False)
    if sample_exists:
      print("Success: Overwritting %s" %sample_file)
    else:
      print("Success: Written to %s" %sample_file)
  else:
    print("OK: File allready exists, skipping! Disable the keep_configs ´-C´ to avoid this.")
    return(False)

  return(True)


def write_subsample_sheet(subsample_sheet, pepdir, keep):
  subsample_file = f"{pepdir}/subsample_sheet.csv"
  subsample_exists = os.path.exists(subsample_file)

  print("Writting subsample sheet", end = "... ")
  if not subsample_exists or not keep:
    subsample_sheet.to_csv(subsample_file, index = False)

    if subsample_exists:
      print("Success: Overwritting %s" %subsample_file)
    else:
      print("Success: Written to %s" %subsample_file)

  else:
    print("OK: File allready exists, skipping! Disable the keep_configs ´-C´ to avoid this.")
    return(False)

  return(True)


def write_PEP(pepdir, keep):
  # Define pep files
  pep_file = f"{pepdir}/project_config.yaml"
  pep_exists = os.path.exists(pep_file)
  
  if not pep_exists or not keep:
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
      print("Success: Written to %s" %pep_file)
    elif pep_exists:
      print("Success: Overwriting %s" %pep_file)
      
  else:
    print("Warning: File allready exists, skipping! Disable the keep_configs ´-C´ to avoid this.")
  

def generate_sheets(reads_dir, assembly_dir, blacklisting, blacklist_file, outdir, tmpdir, keep):
  # Generating dirs
  os.makedirs(outdir, exist_ok=True)
  os.makedirs(tmpdir, exist_ok=True)

  # Screening input files
  reads_metadata = screen_files(directory = reads_dir, type = "Reads")
  assembly_metadata = screen_files(directory = assembly_dir, type = "Assembly")

  metadata = pandas.concat([reads_metadata, assembly_metadata], ignore_index = True)

  # Inspect blacklist if enabled
  blacklist_exists = os.path.exists(blacklist_file)  
  if blacklisting and blacklist_exists:

    blacklist = pandas.read_csv(blacklist_file, sep = "\t")
    exclude_samples = blacklist["file"].values
  
    # Filter metadata
    metadata = metadata[~metadata["file"].isin(exclude_samples)].reset_index()
  
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
    sample_sheet_updated = write_sample_sheet(sample_sheet, pepdir, keep = keep)

    subsample_sheet = metadata[["sample_name", "mate", "file"]]
    subsample_sheet_updated = write_subsample_sheet(subsample_sheet, pepdir, keep = keep)
    sample_files = metadata['file']

    # Generate PEP configuration files.
    write_PEP(pepdir, keep)
  else:
    print("No new samples detected after blacklisting.")
    sample_files = []

  return(sample_files)


def update_blacklist(blacklisting, blacklist_file, clean, sample_files):
  blacklist_exists = os.path.isfile(blacklist_file)

  if blacklisting and clean or not blacklist_exists:
    sample_files.to_csv(blacklist_file, sep = "\t", index=False)
    return(True)
  elif blacklisting:
    sample_files.to_csv(blacklist_file, sep = "\t", index=False, mode = "a", header = False)
  
  return(False)


# Derrive arguments
args = parse_arguments()

reads_dir = args.reads_dir
assembly_dir = args.assembly_dir
database = args.database
threshold = args.threshold
outdir = args.outdir
blacklisting = args.blacklisting
clean_blacklisting = args.clean_blacklisting
multithreading = args.multithreading
keep_tmp = args.keep_tmp
keep_configs = args.keep_configs
threads = args.threads
force_results = args.force_results
force = args.force
skipmake = args.skipmake
dry_run = args.dry_run
debug = args.debug
tmpdir = f"{outdir}/tmp"
blacklist_file = f"{outdir}/blacklist.tsv"

# Validate snakemake structure
validate_snakemake(debug)

# Prepare config file for snakemake
generate_configfile(database, threshold, outdir, tmpdir, multithreading, threads, debug, skipmake)

# Generate subsample sheet
sample_files = generate_sheets(reads_dir, assembly_dir, blacklisting, blacklist_file, outdir, tmpdir, keep = keep_configs)

if len(sample_files) == 0:
  print("Nothing to do exitting!")
  sys.exit(0)

if skipmake:
  print("Setup finished.\nWarning: Skipping Snakemake. To run the pipeline don't use the -S flag!")
else:
  snake_args = ""
  if force:
    snake_args += " -F "
  elif force_results:
    snake_args += " --forcerun all "
  if dry_run:
    snake_args += " -n "

  snakemake_cmd = "snakemake --use-conda --cores %s%s" %(threads, snake_args) 
  if debug:
    print("Running command: %s" %snakemake_cmd)

  subprocess.Popen(snakemake_cmd, shell = True).wait() #check_call() Try except CalledProcessError

  update_blacklist(blacklisting, blacklist_file, clean_blacklisting, sample_files)


if not keep_tmp:
  print("Cleaning up temporary files.")
  shutil.rmtree(f"{outdir}/tmp")

print("All Done!")

