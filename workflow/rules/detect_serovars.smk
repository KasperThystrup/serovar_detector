rule detect_capsules:
	input:
		mate1 = "%s/Read_links/{sample}_R1.fastq.gz" %outdir,
		mate2 = "%s/Read_links/{sample}_R2.fastq.gz" %outdir
	params:
		db = database,
		prefix  = "%s/kma/{sample}/{sample}" %outdir
	output:
		kma_dir = directory("%s/kma/{sample}/" %outdir),
		res_file = "%s/kma/{sample}/{sample}.res" %outdir
	conda:
		"../envs/kma.yaml"
	threads:
		#workflow.cores #kma_threads
		1
	shell:
		"""
		mkdir -p {output.kma_dir}
		#kma -ipe {input.mate1} {input.mate2} -o {params.prefix} -t_db {params.db} -1t1 -dense -ref_fsa -cge -ef -t {threads}
		kma -ipe {input.mate1} {input.mate2} -o {params.prefix} -t_db {params.db} -t {threads}
		"""
