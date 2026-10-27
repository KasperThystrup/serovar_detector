rule copy_database:
  input:
  	database = database
  output:
    fasta = f"{database_prefix}.fasta"
  message:
  	f"Copying {database_name} to system database folder"
  shell:
    """
    cp {input.database} {output.fasta}
    """


rule generate_kma_index:
  input:
    fasta = rules.copy_database.output.fasta
  params:
    prefix = database_prefix
  output:
  	name = f"{database_prefix}.name"
  conda:
    "../envs/serovar_dependencies.yaml"
  message:
  	f"Creating KMA index for {database_name}"
  shell:
  	"""
  	kma index -i {input.fasta} -o {params.prefix}
  	"""