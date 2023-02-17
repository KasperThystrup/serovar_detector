rule summarize_serovars:
	input:
		kma_dir = expand(rules.detect_capsules.output.kma_dir, sample = samples)
	params:
		threshold_id = threshold_id,
		threshold_cov = threshold_cov,
		debug = debug
	output:
		serovar_file = "%s/Results/serovar.tsv" %outdir
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/summarize_serovars.R"
