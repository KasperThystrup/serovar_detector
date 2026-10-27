rule detect_assembly_capsules:
  input:
    database = rules.copy_database.output.fasta,
    assembly = f"{tmpdir}/assemblies/{{sample}}.fasta"
  output:
    tsv = f"{tmpdir}/blast/{{sample}}.tsv"
  conda:
    "../envs/blast.yaml"
  threads:
    1
  message:
    """
    OUTDIR=$(dirname {output.tsv})
    mkdir -p $OUTDIR

    printf 'Template\tAlignment_count\tAlignment_length\tTemplate_length\n' > {output.tsv}

    blastn -subject {input.database} -query {input.assembly} -outfmt '6 sseqid nident length slen' >> {output.tsv}
    """
  shell:
    """
    OUTDIR=$(dirname {output.tsv})
    mkdir -p $OUTDIR

    printf 'Template\tAlignment_count\tAlignment_length\tTemplate_length\n' > {output.tsv}

    blastn -subject {params.db} -query {input.assembly} -outfmt '6 sseqid nident length slen' >> {output.tsv}
    """


rule detect_reads_capsules:
  input:
    database = rules.generate_kma_index.output.name,
    mate1 = f"{tmpdir}/reads/{{sample}}_R1.fastq.gz",
    mate2 = f"{tmpdir}/reads/{{sample}}_R2.fastq.gz"
  params:
    db_prefix = rules.generate_kma_index.params.prefix,
    kma_dir = f"{tmpdir}/kma",
    res_prefix = f"{tmpdir}/kma/{{sample}}"
  output:
    res = f"{tmpdir}/kma/{{sample}}.res"
  conda:
    "../envs/kma.yaml"
  threads:
    1
  message:
    """
    kma -ipe {input.mate1} {input.mate2} -o {params.res_prefix} -t_db {params.db_prefix} -t {threads}
    """
  shell:
    """
    kma -ipe {input.mate1} {input.mate2} -o {params.res_prefix} -t_db {params.db_prefix} -t {threads}
    """

