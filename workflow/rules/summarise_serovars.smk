rule summarize_serovars:
	input:
		assembly_results_dir = expand(rules.detect_assembly_capsules.output.kma_dir, sample = pep.sample_table["sample_name"]),
		assembly_results = expand(rules.detect_assembly_capsules.output.res_file, sample = pep.sample_table["sample_name"]),
		reads_results_dir = expand(rules.detect_reads_capsules.output.kma_dir, sample = pep.sample_table["sample_name"]),
		reads_results = expand(rules.detect_reads_capsules.output.res_file, sample = pep.sample_table["sample_name"])
	params:
		threshold = threshold,
		blacklisting = blacklisting,
		debug = debug
	output:
		serovar_file = "%s/serovar.tsv" %outdir
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/summarize_serovars.R"
