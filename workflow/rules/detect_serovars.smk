rule detect_assembly_capsules:
	input:
		assemblies = "%s/{sample}.fasta" %assembly_dir
	params:
		db = database,
		prefix  = "%s/kma/{sample}/{sample}" %outdir
	output:
		kma_dir = temp(directory("%s/kma/{sample}/" %outdir)),
		res_file = "%s/kma/{sample}/{sample}.res" %outdir
	conda:
		"../envs/kma.yaml"
	threads:
		1
	message:
		"""
                mkdir -p {output.kma_dir}
                kma -i {input.assemblies} -o {params.prefix} -t_db {params.db} -t {threads}
                """
	shell:
		"""
		mkdir -p {output.kma_dir}
		kma -i {input.assemblies} -o {params.prefix} -t_db {params.db} -t {threads}
		"""

rule detect_reads_capsules:
	input:
		mate1 = "%s/{sample}_R1.fastq.gz" %reads_dir,
		mate2 = "%s/{sample}_R2.fastq.gz" %reads_dir
	params:
		db = database,
		prefix  = "%s/kma/{sample}/{sample}" %outdir
	output:
		kma_dir = temp(directory("%s/kma/{sample}/" %outdir)),
		res_file = "%s/kma/{sample}/{sample}.res" %outdir
	conda:
		"../envs/kma.yaml"
	threads:
		1
	message:
		"""
                mkdir -p {output.kma_dir}
                kma -ipe {input.mate1} {input.mate2} -o {params.prefix} -t_db {params.db} -t {threads}
                """
	shell:
		"""
		mkdir -p {output.kma_dir}
		kma -ipe {input.mate1} {input.mate2} -o {params.prefix} -t_db {params.db} -t {threads}
		"""

