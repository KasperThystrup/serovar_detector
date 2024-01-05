rule link_assemblies:
    input:
        assembly=lambda wildcards: assembly_sheet.loc[wildcards.sample, 'file']
    output:
        assembly="%s/tmp/{sample}.fasta" % outdir
    shell:
        """
        ln -s '{input.assembly}' '{output.assembly}'
        """


rule link_reads:
    input:
        mate1=lambda wildcards: reads_sheet.loc[wildcards.sample, 'file'][0],
        mate2=lambda wildcards: reads_sheet.loc[wildcards.sample, 'file'][1]
    output:
        mate1="%s/tmp/{sample}_R1.fastq.gz" % outdir,
        mate2="%s/tmp/{sample}_R2.fastq.gz" % outdir
    shell:
        """
        ln -s '{input.mate1}' '{output.mate1}'
        ln -s '{input.mate2}' '{output.mate2}'
        """

