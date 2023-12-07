rule detect_assembly_capsules:
	input:
		assemblies = lambda assembly_path: pep.sample_table.loc[assembly_path.sample]["assembly"]
	params:
		db = database,
		prefix  = lambda kma_path: "%s/%s" %(pep.sample_table.loc[kma_path.sample]["kma"], os.path.basename(pep.sample_table.loc[kma_path.sample]["kma"]))
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
		mate1 = lambda read1_path: pep.sample_table.loc[read1_path.sample]["read1"],
                mate2 = lambda read2_path: pep.sample_table.loc[read2_path.sample]["read2"]
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

