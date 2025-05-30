rule summarize_serovars:
	input:
		assembly_results = expand(rules.detect_assembly_capsules.output.res_file, sample = assembly_sheet["sample_name"].values.tolist()),
		reads_results = expand(rules.detect_reads_capsules.output.res_file, sample = reads_sheet["sample_name"].values.tolist())
	params:
		threshold = threshold,
		debug = debug
	output:
		serovar_file = "%s/serovar.tsv" %outdir
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/summarize_serovars.R"

