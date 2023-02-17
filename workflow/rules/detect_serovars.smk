rule detect_capsules:
	input:
		mate1 = "%s/{sample}_R1.fastq.gz" %sample_dir,
		mate2 = "%s/{sample}_R2.fastq.gz" %sample_dir
	params:
		db = database,
		prefix  = "%s/kma/{sample}/{sample}" %outdir
	output:
		kma_dir = temp(directory("%s/kma/{sample}/" %outdir)),
		res_file = "%s/kma/{sample}/{sample}.res" %outdir
	conda:
		"../envs/kma.yaml"
	threads:
		kma_threads
	shell:
		"""
		mkdir -p {output.kma_dir}
		#kma -ipe {input.mate1} {input.mate2} -o {params.prefix} -t_db {params.db} -1t1 -dense -ref_fsa -cge -ef -t {threads}
		kma -ipe {input.mate1} {input.mate2} -o {params.prefix} -t_db {params.db} -t {threads}
		"""
