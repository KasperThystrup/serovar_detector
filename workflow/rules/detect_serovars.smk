rule detect_assembly_capsules:
  input:
    assembly_file = rules.link_assemblies.output.assembly
  params:
    db = database,
    prefix = lambda wildcards: assembly_seet.loc[wildcards.sample, "sample"]
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
    kma -i {input.assembly_file} -o {params.prefix} -t_db {params.db} -t {threads}
    """
  shell:
    """
    mkdir -p {output.kma_dir}
    kma -i {input.assembly_file} -o {params.prefix} -t_db {params.db} -t {threads}
    """


rule detect_reads_capsules:
  input:
    mate1 = rules.link_reads.output.mate1,
    mate2 = rules.link_reads.output.mate2
  params:
    db = database,
    prefix = lambda wildcards: reads_sheet.loc[wildcards.sample, "sample"]
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

