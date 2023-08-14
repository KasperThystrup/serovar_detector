rule summarize_serovars:
	input:
		assembly_results_dir = expand(rules.detect_assembly_capsules.output.kma_dir, sample = sample_assemblies),
		assembly_results = expand(rules.detect_assembly_capsules.output.res_file, sample = sample_assemblies),
		reads_results_dir = expand(rules.detect_reads_capsules.output.kma_dir, sample = sample_reads),
		reads_results = expand(rules.detect_reads_capsules.output.res_file, sample = sample_reads)
	params:
		threshold = threshold,
		debug = debug
	output:
		serovar_file = "%s/serovar.tsv" %outdir
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/summarize_serovars.R"
