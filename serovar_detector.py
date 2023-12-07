import argparse
import os
import sys
import yaml
import glob
import subprocess
import re
import pandas

def parse_arguments():
  parser = argparse.ArgumentParser(description = "Screen read files and assemblies for Serovar biomarker genes, in order to preovide suggestions for isolate serovar. Currently only supporting Actinobacillus Pleuropneumoniae.")
  parser.add_argument("-r", metavar = "--reads_dir", dest = "reads_dir", help = "Input path to reads directory", required = False)
  parser.add_argument("-a", metavar = "--assembly_dir", dest = "assembly_dir", help = "Input path to assembly directory", required = False)
  parser.add_argument("-D", metavar = "--database", dest = "database", help = "Path and prefix to kmer-aligner database", required = True)
  parser.add_argument("-T", metavar = "--theshold", dest = "threshold", help = "Cutoff threshold of match coverage and identity. Ignore threshold by setting to 0 or False. (Default 98)", default = 98)
  parser.add_argument("-o", metavar = "--outdir", dest = "outdir", help = "Output path to Results and Temporary files directory", required = True)
  parser.add_argument("-b", dest = "blacklisting", help = "Blacklist succesfsfully analysed samples, usable for surveillance / continous projects. (Default False)", action = "store_true")
  parser.add_argument("-B", dest = "clean_blacklist", help = "Remove existing blacklist file. (Default False)", action = "store_true")
  parser.add_argument("-M", dest = "multithreading", help = "Enable multithreading during kmer alignment, use for huge samples only. (Default False)", action = "store_true")
  parser.add_argument("-k", dest = "keep", help = "Preserve temporary files such as KMA result files. (Default False)", action = "store_true")
  parser.add_argument("-t", metavar = "--threads", dest = "threads", help = "Number of threads to allocate for the pipeline. (Default 1)", default = 1)
  parser.add_argument("-f", dest = "force_results", help = "Force rerun of the Results tasks in the pipeline. (Default False)", action = "store_true")
  parser.add_argument("-F", dest = "force", help = "Force rerun of all tasks in pipeline. (Default False)", action = "store_true")
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


def generate_configfile(reads_dir, assembly_dir, database, threshold, outdir, blacklisting, multithreading, threads, debug, skipmake):
  # Define config file
  config_file = "config/config.yaml"
  
  # Check for existing config file
  if not os.path.isfile(config_file):
    config_dir = os.path.dirname(config_file)
    
    # Check for existing config dir
    if not os.path.isdir(config_dir):
      print("No config dir detected, creating directory.")
      os.mkdir(config_dir)

  # Directing full paths
  if reads_dir is not None:
    reads_dir_path = os.path.abspath(reads_dir).rstrip("/")
  else:
    reads_dir_path = ""
  if assembly_dir is not None:
    assembly_dir_path = os.path.abspath(assembly_dir).rstrip("/")
  else:
    assembly_dir_path = ""

  out_path = os.path.abspath(outdir).rstrip("/")
  
  config = {
    "database": database, "threshold": threshold,  "outdir": out_path,
    "blacklisting": blacklisting, "multithreading": multithreading, "debug": debug
  }

  with open(config_file, "w") as config_yaml:
    yaml.dump(config, config_yaml)


def generate_subsample_sheet(reads_dir, assembly_dir):
  
  print("Screening assembly directory for files", end = "... ")
  assembly_sheet = pandas.DataFrame(columns = ["sample_name", "file", "ext"])

  if assembly_dir is not None and os.path.exists(assembly_dir):
    # Screen for files
    sample_fasta = [sample for sample in glob.glob("%s/*fasta" %assembly_dir)]
    
    # Generate grouped search object
    assembly_search = [
      re.search(
        "^(?P<file>\S+\/(?P<sample_name>\S+)\.(?P<ext>[fast]+))" , fasta_file
      ) for fasta_file in sample_fasta
    ]
    
    # Generate DataFrame from search object groups
    assembly_raw = pandas.DataFrame({
      "sample_name": [search.group("sample_name") for search in assembly_search],
      "file": [search.group("file") for search in assembly_search],
      "ext": [search.group("ext") for search in assembly_search]
    })
    assembly_sheet = assembly_raw.sort_values(by = "sample_name")

    print("Success: %s assembly files found" %len(sample_fasta))
  # Report on screening status
  elif assembly_dir is None:
    print("OK: No assembly directory provided, skipping!")
  else:
    print("Failed: Assembly directory does not exist!\n  - %s" %assembly_dir)

  print("Screening reads directory for files", end = "... ")
  reads_sheet = pandas.DataFrame(columns = ["sample_name", "mate", "file", "ext"])
  
  if reads_dir is not None and os.path.exists(reads_dir):
    # Screen for files
    sample_fastq = [sample_read for sample_read in glob.glob("%s/*fastq*" %reads_dir)]
      
    # Generate grouped search object
    reads_search = [
      re.search(
        "^(?P<file>\S+/(?P<sample_name>\S+?)(?:_S\d+)?(?:_L\d+)?_(?P<mate>R[12])(?:_\d+)?\.(?P<ext>fastq\.gz))",
        fastq_file
      ) for fastq_file in sample_fastq
    ] 
    
    # Generate DataFrame from search object groups
    reads_raw = pandas.DataFrame({
      "sample_name": [search.group("sample_name") for search in reads_search],
      "mate": [search.group("mate") for search in reads_search],
      "file": [search.group("file") for search in reads_search],
      "ext": [search.group("ext") for search in reads_search]
    })
    reads_sheet = reads_raw.sort_values(by = ["sample_name", "mate"])

    print("Success: %s read files found" %len(sample_fastq))
  # Report on screening status
  elif reads_dir is None:
    print("OK: No reads directory provided, skipping!")
  else:
    print("Failed: Reads directory does not exist!\n  - %s" %reads_dir)

  subsample_concattenated = pandas.concat([reads_sheet, assembly_sheet], sort = False)
  subsample_sorted = subsample_concattenated.sort_values(by = ["ext", "sample_name"], ascending=False)

  return(subsample_sorted[["sample_name", "mate", "file", "ext"]])


def generate_sample_sheet(subsample_sheet):
  print("Generating sample sheet", end = "... ")

  file_count = len(subsample_sheet.index)
  
  if file_count > 0:
    sample_subset = subsample_sheet[["sample_name"]]
    sample_sheet = sample_subset.drop_duplicates()

    print("Success: A total of %s samples have been annotated!" %len(sample_sheet.index))
    return(sample_sheet)
  else:
    print("Failed: Subsample sheet has no rows")
    return(False)


def write_subsample_sheet(subsample_sheet, force):
  subsample_file = "schemas/subsample_sheet.csv"
  subsample_exists = os.path.exists(subsample_file)

  print("Writting subsample sheet", end = "... ")
  if not subsample_exists or force:
    subsample_sheet.to_csv(subsample_file, index = False)

    if not force:
      print("Success: Written to %s" %subsample_file)
    else:
      print("Success: Overwritting %s" %subsample_file)

  else:
    print("OK: File allready exists, skipping! To renew, delete/rename the old file or enable the `-f` (force) option.")
    return(False)

  return(True)


def write_sample_sheet(sample_sheet, force):
  sample_file = "schemas/sample_sheet.csv"
  sample_exists = os.path.exists(sample_file)

  print("Writting sample sheet", end = "... ")
  if not sample_exists or force:
    sample_sheet.to_csv(sample_file, index = False)
    if not force:
      print("Success: Written to %s" %sample_file)
    else:
      print("Success: Overwirtting %s" %sample_file)
  else:
    print("OK: File allready exists, skipping! To renew, delete/rename the old file or enable the `-f` (force) option.")
    return(False)

  return(True)


def write_PEP(subsample_updated, outdir, force):
   # Generate PEP configuration:
  PEP_header = "pep_version: 2.1.0\n"
  PEP_sample = "sample_table: 'sample_sheet.csv'\n"
  PEP_subsample = "subsample_table: 'subsample_sheet.csv'"
  PEP_modifiers = """
sample_modifiers:
  imply:
    - if:
        ext: ["fastq.gz", "fastq.gz"]
      then:
        type: "reads"
    - if:
        ext: "fasta"
      then:
        type: "assembly"
"""
  
  pep_file = "schemas/project_config.yaml"
  pep_exists = os.path.exists(pep_file)
  print("Writing project configuration file", end = "... ")
  if not pep_exists or subsample_updated or force:
    
    with open(pep_file, "w") as config_file:
      config_file.write(PEP_header)
      config_file.write(PEP_sample)
      config_file.write(PEP_subsample)
      config_file.write(PEP_modifiers)
    
    if not pep_exists:
      print("Success: Written to %s" %pep_file)
    elif pep_exists and subsample_updated and not force:
      print("Success: Subsample sheet updated! Overwritting %s" %pep_file)
    else:
      print("Success: Overwriting %s" %pep_file)
  else:
    print("OK: File allready exists, skipping! To renew, delete/rename the old file or enable the `-f` (force) option.")
    return(False)
  
  return(True)

 

# Derrive arguments
args = parse_arguments()

reads_dir = args.reads_dir
assembly_dir = args.assembly_dir
database = args.database
threshold = args.threshold
outdir = args.outdir
blacklisting = args.blacklisting
clean_blacklist = args.clean_blacklist
multithreading = args.multithreading
keep = args.keep
threads = args.threads
force_results = args.force_results
force = args.force
skipmake = args.skipmake
dry_run = args.dry_run
debug = args.debug

# Validate snakemake structure
validate_snakemake(debug)

# Prepare config file for snakemake
generate_configfile(reads_dir, assembly_dir, database, threshold, outdir, blacklisting, multithreading, threads, debug, skipmake)


# Preparing output directory
if not os.path.exists("schemas"):
  print("Creating PEP directory:")
  os.mkdir("schemas")

# Generate subsample sheet
subsample_sheet = generate_subsample_sheet(reads_dir, assembly_dir)
sample_sheet = generate_sample_sheet(subsample_sheet)
subsample_updated = write_subsample_sheet(subsample_sheet, force)
sample_updated = write_sample_sheet(sample_sheet, force)

pep_updated = write_PEP(subsample_updated, outdir, force)
sys.exit(0)
if skipmake:
  print("Warning: Skipping Snakemake!")
else:
  snake_args = ""
  if keep:
    snake_args += " --notemp "
  if force:
    snake_args += " -F "
  elif force_results:
    snake_args += " --forcerun all "
  if dry_run:
    snake_args += " -n "
  if clean_blacklist:
    snake_args += "--forcerun clean_blacklist"
  

  snakemake_cmd = "snakemake --use-conda --cores %s%s" %(threads, snake_args) 
  if debug:
    print("Running command: %s" %snakemake_cmd)

  subprocess.Popen(snakemake_cmd, shell = True).wait()

