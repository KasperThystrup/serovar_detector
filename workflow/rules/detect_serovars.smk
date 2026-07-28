rule detect_assembly_capsules:
  input:
    assembly = "%s/assemblies/{sample}.fasta" %tmpdir
  params:
    db = blast_database
  output:
    tsv = "%s/blast/{sample}.tsv" %tmpdir
  conda:
    "../envs/blast.yaml"
  threads:
    1
  message:
    """
    OUTDIR=$(dirname {output.tsv})
    mkdir -p $OUTDIR

    printf 'Template\tAlignment_count\tAlignment_length\tTemplate_length\n' > {output.tsv}

    blastn -subject {params.db} -query {input.assembly} -outfmt '6 sseqid nident length slen' >> {output.tsv}
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
    mate1 = "%s/reads/{sample}_R1.fastq.gz" %tmpdir,
    mate2 = "%s/reads/{sample}_R2.fastq.gz" %tmpdir
  params:
    db = database,
    prefix = "%s/kma/{sample}" %tmpdir,
    kma_dir = "%s/kma" %tmpdir
  output:
    res = "%s/kma/{sample}.res" %tmpdir
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

