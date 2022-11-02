rule serovar_profiles:
	params:
		debug = debug
	output:
		profiles_file = "%s/profiles/serovar_profiles.tsv" %outdir
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/generate_serovar_profiles.R"


rule summarize_serovars:
	input:
		kma_dir = expand(rules.detect_capsules.output.kma_dir, sample = samples),
		profiles_file = rules.serovar_profiles.output.profiles_file
	params:
		threshold_id = threshold_id,
		threshold_cov = threshold_cov,
		debug = debug
	output:
		serovar_file = "%s/serovars.tsv" %outdir
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/summarize_serovars.R"