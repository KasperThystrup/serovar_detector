rule summarize_serovars:
  input:
    assembly_files = expand(rules.detect_assembly_capsules.output.res_file, sample = assembly_sheet["sample_name"].values.tolist()),
    reads_files = expand(rules.detect_reads_capsules.output.res_file, sample = reads_sheet["sample_name"].values.tolist())
  params:
    threshold = threshold,
    serovar_profiles = serovar_profiles,
    debug = debug
  output:
    results_wide_file = "%s/serovar.tsv" %outdir,
    results_long_file = "%s/results_long.tsv" %outdir
  conda:
    "../envs/R.yaml"
  script:
    "../scripts/summarize_serovars.R"

