import argparse
import os
import sys
import yaml
import glob
import subprocess

def parse_arguments():
  parser = argparse.ArgumentParser(description = "Screen assemblies for junction sequences in order to characterize specific plasmid markers")
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
  
  config = {"reads_dir" : reads_dir_path, "assembly_dir" : assembly_dir_path, "database" : database, "threshold" : threshold,  "outdir" : out_path, "blacklisting" : blacklisting, "multithreading" : multithreading, "debug" : debug}

  with open(config_file, "w") as config_yaml:
    yaml.dump(config, config_yaml)

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
  if blacklisting:
    snake_args += " --forcerun blacklist "
  if clean_blacklist:
    snake_args += "--forcerun clean_blacklist"
  

  snakemake_cmd = "snakemake --use-conda --cores %s%s" %(threads, snake_args) 
  if debug:
    print("Running command: %s" %snakemake_cmd)

  subprocess.Popen(snakemake_cmd, shell = True).wait()

