rule detect_assembly_capsules:
  input:
    assembly_file = "%s/assemblies/{sample}.fasta" %outdir
  params:
    db = database,
    prefix =  "%s/kma/{sample}" %outdir,
    kma_dir = directory("%s/kma" %outdir)
  output:
    res_file = "%s/kma/{sample}.res" %outdir
  conda:
    "../envs/kma.yaml"
  threads:
    1
  message:
    """
    mkdir -p {params.kma_dir}
    kma -i {input.assembly_file} -o {params.prefix} -t_db {params.db} -t {threads}
    """
  shell:
    """
    mkdir -p {params.kma_dir}
    kma -i {input.assembly_file} -o {params.prefix} -t_db {params.db} -t {threads}
    """


rule detect_reads_capsules:
  input:
    mate1 = "%s/reads/{sample}_R1.fastq.gz" %outdir,
    mate2 = "%s/reads/{sample}_R2.fastq.gz" %outdir
  params:
    db = database,
    prefix = "%s/kma/{sample}" %outdir,
    kma_dir = directory("%s/kma" %outdir)
  output:
    res_file = "%s/kma/{sample}.res" %outdir
  conda:
    "../envs/kma.yaml"
  threads:
    1
  message:
    """
    mkdir -p {params.kma_dir}
    kma -ipe {input.mate1} {input.mate2} -o {params.prefix} -t_db {params.db} -t {threads}
    """
  shell:
    """
    mkdir -p {params.kma_dir}
    kma -ipe {input.mate1} {input.mate2} -o {params.prefix} -t_db {params.db} -t {threads}
    """

