rule detect_assembly_capsules:
  input:
    assembly_file = "%s/assemblies/{sample}.fasta" %tmpdir
  params:
    db = database,
    prefix =  "%s/assemblies/kma/{sample}" %tmpdir,
    kma_dir = "%s/assemblies/kma" %tmpdir
  output:
    res_file = "%s/assemblies/kma/{sample}.res" %tmpdir
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
    mate1 = "%s/reads/{sample}_R1.fastq.gz" %tmpdir,
    mate2 = "%s/reads/{sample}_R2.fastq.gz" %tmpdir
  params:
    db = database,
    prefix = "%s/reads/kma/{sample}" %tmpdir,
    kma_dir = "%s/reads/kma" %tmpdir
  output:
    res_file = "%s/reads/kma/{sample}.res" %tmpdir
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

