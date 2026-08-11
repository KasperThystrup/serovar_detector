read_kma <- function(res_file){
  
  logger::log_debug("Performing pattern matching to determine sample name from file: ", res_file)
  file_matches <- res_file |>
    basename() |>
    stringr::str_match_all(pattern = "^(?<sample>[^\\.]+(\\.[^\\.]+)*)\\.(?<ext>\\w+)$") |>
    as.data.frame()

  sample_name <- file_matches$sample
  
  logger::log_debug("Reading and updating headers.")
  header <- res_file |>
    readr::read_lines(n_max = 1) |>
    strsplit(split = "\t") |>
    unlist()
  
  header[1] <- "Template"

  logger::log_debug("Reading file.\n\n")
  readr::read_tsv(file = res_file, col_names = header, col_types = "ciiiddddddd", skip = 1) |>
    tidyr::separate(
      col = Template,
      into = c("Template_Gene", "Template_Strain", "Template_Serovar"),
      sep = "_"
    ) |>
    # Adding sample name column
    dplyr::mutate(Sample = sample_name, Mapper = "KMA") |>
    dplyr::select(Sample, Mapper, Template_Gene, Template_Strain, Template_Serovar, Template_Identity, Template_Coverage)
}


read_blast <- function(tsv_file){
  logger::log_debug("Performing pattern matching to determine sample name from file: ", tsv_file)
  file_matches <- tsv_file |>
    basename() |>
    stringr::str_match_all(pattern = "^(?<sample>[^\\.]+(\\.[^\\.]+)*)\\.(?<ext>\\w+)$") |>
    as.data.frame()
  
  sample_name <- file_matches$sample
  
  logger::log_debug("Reading file.\n\n")
  blast_raw <- readr::read_tsv(file = tsv_file, col_types = "ciii") |>
    dplyr::mutate(
      Sample = sample_name,
      Mapper = "Blastn",
      .before = 1
    )
  
  if (nrow(blast_raw) == 0){
    blast <- tibble::tibble(
      Sample = character(),
      Mapper = character(),
      Template_Gene = character(),
      Template_Strain = character(),
      Template_Serovar = character(),
      Template_Identity = numeric(),
      Template_Coverage = numeric()
    )
    
  } else {
    blast_details <- tidyr::separate(
      blast_raw,
      col = Template,
      into = c("Template_Gene", "Template_Strain", "Template_Serovar"),
      sep = "_"
    ) |>
      dplyr::group_by(Sample, Template_Gene) |>
      dplyr::filter(Alignment_length == max(Alignment_length)) |>
      dplyr::filter(Alignment_length == max(Alignment_length)) |>
      dplyr::slice(1)

    # Adding sample name column
    blast <- dplyr::mutate(
      blast_details,
      Template_Coverage = Alignment_length / Template_length * 100,
      Template_Identity = Alignment_count / Template_length * 100
    ) |>
      dplyr::select(Sample, Mapper, Template_Gene, Template_Strain, Template_Serovar, Template_Identity, Template_Coverage)
    
  }
  
  return(blast)
}



apply_thresholds <- function(res_table, threshold){
  logger::log_debug("Applying thresholds")
  
  empty_table <- nrow(res_table) == 0
  
  if (empty_table)
    return(res_table)
  
  dplyr::mutate(
    res_table,
    match_perfect = Template_Identity == 100 & Template_Coverage == 100,
    match_imperfect = !match_perfect & (Template_Identity >= threshold & Template_Coverage >= threshold), # Manual thresholds!!
    match_partial = !match_imperfect & !match_perfect
  )
}


object_to_dataframe <- function(list_iterator, object){
  item <- object[list_iterator]
  serovar <- names(item)
  genes <- unlist(item)
  
  dplyr::tibble(Gene = genes, Serovar = serovar)
}


list_to_table <- function(object){
  purrr::map_dfr(seq_along(object), object_to_dataframe, object)
}


generate_serovar_profiles <- function(serovar_config_yaml){
  logger::log_debug("Importing serovar profiles from: ", serovar_config_yaml)
  serovar_yaml <- yaml::yaml.load_file(input = serovar_config_yaml)
  
  logger::log_debug("Converting format into dataframe.")
  list_to_table(serovar_yaml)
}


resolve_serovars <- function(results, serovar_profiles){
  
  logger::log_debug("Merging kma results with serovar profiles table.") ### Warnings!!!
  profiles <- dplyr::select(
    results, 
    Sample, Mapper, Template_Gene, Template_Strain, Template_Serovar, Template_Identity,
    Template_Coverage, match_perfect, match_imperfect, match_partial
  ) |>
    dplyr::left_join(
      y = serovar_profiles, by = c("Template_Gene" = "Gene"), relationship = "many-to-many"
    ) |>
    dplyr::group_by(Sample, Mapper, Serovar)
  
  logger::log_debug("Counting the capsule genes to determine repressentation of serovars.")
  overview <- dplyr::reframe(
    profiles,
    Gene_count = sum(c(match_perfect, match_imperfect))
  ) |>
    
    # If Mapper is removed here, it will be affected by winner-takes-it-all
    dplyr::group_by(Sample, Mapper) |>
    # Determine which serovars are best represented relative to capsule gene counts
    dplyr::reframe(
      Serovar,
      Gene_count,
      winner = Gene_count == max(Gene_count)
    )
  
  logger::log_debug("Counting expected amount of capsule genes for each serovar")
  profiles_count <- dplyr::group_by(serovar_profiles, Serovar) |>
    dplyr::summarise(capsule_count = dplyr::n(), .groups = "keep") ## .groups added
  
  logger::log_debug("Filtering the most repressented serovar and quantifying capsule gene frequency.")
  serovar_suggestions <- dplyr::filter(overview, winner) |>
    dplyr::group_by(Sample, Mapper) |>
    dplyr::summarise(
      suggestions = dplyr::n(),
      Serovar = paste(Serovar, collapse = ","),
      count = dplyr::case_when(
        suggestions != 1 ~ NA_integer_,
        TRUE ~ unique(Gene_count)
      ),
      .groups = "keep"
    )
  
  serovars <- dplyr::left_join(
    x = serovar_suggestions, y = profiles_count, by = "Serovar"
  ) |>
    dplyr::summarise(
      Sample,
      Suggested_serovar = dplyr::case_when(
        is.na(count) ~ paste(Serovar, sep = ", "),
        count / capsule_count < 0.5 ~ paste0("?", Serovar, collapse = ", "),
        TRUE ~ Serovar
      ),
        
      Frequency = dplyr::case_when(
        is.na(count) ~ NA_character_,
        TRUE ~ paste(count, capsule_count, sep = " of ")
      ),
      .groups = "keep"
    )

  merged <- dplyr::left_join(profiles, serovars, by = c("Sample", "Mapper")) |>
    dplyr::group_by(Sample, Template_Gene)
  
  detailed <- dplyr::mutate(
    merged,
    member = any(Serovar == Suggested_serovar),
    gene_id = dplyr::case_when(
      Template_Identity == 100 ~ "",
      Template_Identity != 100 & Template_Coverage == 100 ~ paste0(Template_Identity, "% ID"),
      TRUE ~ paste0(Template_Identity, "% ID, ")
    ),
    gene_cov = dplyr::case_when(
      Template_Coverage == 100 ~ "",
      TRUE ~ paste0(Template_Coverage, "% COV")
    ),
    gene_detailed = dplyr::case_when(
      match_perfect ~ Template_Gene,
      TRUE ~ paste0(Template_Gene, " (", gene_id, gene_cov, ")")
    ),
    class = stringr::str_extract(string = Template_Gene, pattern = "\\w$")
  ) |>
    dplyr::group_by(Sample, Mapper) |>
    dplyr::arrange(class)
  
  logger::log_debug("Defining and annotating accepted gene- and partial gene-matches.")
  genes_raw <- dplyr::reframe(
    detailed,
    Serovar_match = dplyr::case_when(
      member & match_perfect ~ gene_detailed,
      member & match_imperfect ~ gene_detailed
    ),
    Serovar_partial = dplyr::case_when(
      member & match_partial ~ gene_detailed
    ),
    
    Others_match = dplyr::case_when(
      !member & match_perfect ~ gene_detailed,
      !member & match_imperfect ~ gene_detailed
    ),
    Others_partial = dplyr::case_when(
      !member & match_partial ~ gene_detailed
    )
  ) |>
    dplyr::group_by(Sample, Mapper)
  
  logger::log_debug("Removing NA's and compressing accepted- and partial -genes into single columns.")
  genes <- dplyr::summarise(
    genes_raw,
    Serovar_match = paste0(unique(Serovar_match[!is.na(Serovar_match)]), collapse = ", "),
    Serovar_partial = paste0(unique(Serovar_partial[!is.na(Serovar_partial)]), collapse = ", "),
    Others_match = paste0(unique(Others_match[!is.na(Others_match)]), collapse = ", "),
    Others_partial = paste0(unique(Others_partial[!is.na(Others_partial)]), collapse = ", "),
    .groups = "drop"
  )
  
  dplyr::left_join(serovars, genes, by = c("Sample", "Mapper"))
    
}
  

summarize_serovars <- function(blast_results, kma_results, serovar_config_yaml, threshold, serovar_file){
  
  logger::log_info("Reading and merging all .res files.")
  kma_table <- purrr::map_dfr(kma_results, read_kma)

  blast_table <- purrr::map_dfr(blast_results, read_blast)
  
  all_results <- rbind(kma_table, blast_table)
  
  filtered <- apply_thresholds(all_results, threshold)
  
  if (nrow(filtered) < 1){
    message("Serovar results file(s) are empty!")
    readr::write_file(x = "No serovars", file = serovar_file)

  } else {
    
    logger::log_info("Generating serovar profiles from profile-config file.")
    serovar_profiles <- generate_serovar_profiles(serovar_config_yaml)
    
    logger::log_info("Determining the most frequently repressented serovars and serovar genes.")
    results <- resolve_serovars(filtered, serovar_profiles)
    
    logger::log_info("Writing results to: ", serovar_file)
    readr::write_tsv(x = dplyr::arrange(results, Sample), file = serovar_file)
  }  
}

blast_results <- snakemake@input$blast_results
kma_results <- snakemake@input$kma_results
threshold <- snakemake@params$threshold
serovar_profiles <- snakemake@params$serovar_profiles
dbg <- snakemake@params$debug
serovar_file <- snakemake@output$serovar_file

logger::log_threshold(level = logger::INFO)
if (dbg){
  logger::log_threshold(level = logger::DEBUG)
  
  tmp_file <- file.path(
    dirname(serovar_file),
    "summarize_serovars.RData"
  )
  
  logger::log_debug("Writing snakemake R object to ", tmp_file)
  save(snakemake, file = tmp_file)
  
}

summarize_serovars(
  blast_results = blast_results,
  kma_results = kma_results,
  serovar_config_yaml = serovar_profiles,
  threshold = threshold,
  serovar_file = serovar_file
)

