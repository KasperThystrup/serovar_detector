rule detect_assembly_capsules:
  input:
    assembly_file = "%s/{sample}/{sample}.fasta" %outdir
  params:
    db = database,
    prefix =  "%s/{sample}/kma/{sample}" %outdir
  output:
    kma_dir = directory("%s/{sample}/kma" %outdir),
    res_file = "%s/{sample}/kma/{sample}.res" %outdir
  conda:
    "../envs/kma.yaml"
  threads:
    1
  message:
    """
    mkdir -p {output.kma_dir}
    kma -i {input.assembly_file} -o {params.prefix} -t_db {params.db} -t {threads}
    """
  shell:
    """
    mkdir -p {output.kma_dir}
    kma -i {input.assembly_file} -o {params.prefix} -t_db {params.db} -t {threads}
    """


rule detect_reads_capsules:
  input:
    mate1 = "%s/{sample}/{sample}_R1.fastq.gz" %outdir,
    mate2 = "%s/{sample}/{sample}_R2.fastq.gz" %outdir
  params:
    db = database,
    prefix = "%s/{sample}/kma/{sample}" %outdir
  output:
    kma_dir = directory("%s/{sample}/kma" %outdir),
    res_file = "%s/{sample}/kma/{sample}.res" %outdir
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

