rule summarize_serovars:
	input:
		blast_results = expand(rules.detect_assembly_capsules.output.tsv, sample = assembly_sheet["sample_name"].values.tolist()),
		kma_results = expand(rules.detect_reads_capsules.output.res, sample = reads_sheet["sample_name"].values.tolist())
	params:
		threshold = threshold,
		serovar_profiles = serovar_profiles,
		debug = debug
	output:
		serovar_file = "%s/serovars.tsv" %outdir
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/summarize_serovars.R"

