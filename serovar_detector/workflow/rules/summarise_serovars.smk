rule summarize_serovars:
	input:
		blast_results = expand(rules.detect_assembly_capsules.output.tsv, sample = assembly_samples),
		kma_results = expand(rules.detect_reads_capsules.output.res, sample = reads_samples)
	params:
		threshold = threshold,
		serovar_profiles = serovar_profiles,
		debug = debug
	output:
		serovar_file = f"{outdir}/serovars.tsv"
	conda:
		"../envs/serovar_dependencies.yaml"
	script:
		"../scripts/summarize_serovars.R"

