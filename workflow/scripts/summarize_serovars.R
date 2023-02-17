library(magrittr)


read_res <- function(res_file){
  
  logger::log_debug("Performing pattern matching to determine sample name from file: ", res_file)
  file_matches <- res_file %>%
    basename %>%
    stringr::str_match_all(pattern = "^(?<sample>[^\\.]+)\\.(?<ext>\\w+)$") %>%
    as.data.frame()

  sample_name <- file_matches$sample
  
  logger::log_debug("Reading and updating headers.")
  header <- res_file %>%
    readr::read_lines(n_max = 1) %>%
    strsplit(split = "\t") %>%
    unlist
  
  header[1] <- "Template"
  
  logger::log_debug("Reading file.\n\n")
  readr::read_tsv(file = res_file, col_names = header, col_types = "ciiiddddddd", skip = 1) %>%
    tidyr::separate(
      col = Template,
      into = c("Template_Gene", "Template_Strain", "Template_Serovar"),
      sep = "_"
    ) %>%
    # Adding sample name column
    dplyr::mutate(Sample = sample_name)
}


apply_thresholds <- function(res_table, threshold_id, threshold_cov){
  logger::log_debug("Applying thresholds")
  
  empty_table <- nrow(res_table) == 0
  
  if (empty_table)
    return(NULL)
  
  dplyr::mutate(
    res_table,
    match_perfect = Template_Identity == 100 & Template_Coverage == 100,
    match_imperfect = !match_perfect & (Template_Identity >= threshold_id & Template_Coverage >= threshold_cov), # Manual thresholds!!
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


resolve_serovars <- function(kma_table, profiles){
  
  logger::log_debug("Merging kma results with serovar profiles table.")
  kma_profile <- dplyr::select(
    kma_table, 
    Sample, Template_Gene, Template_Strain, Template_Serovar, Template_Identity,
    Template_Coverage, match_perfect, match_imperfect, match_partial
  ) %>%
    dplyr::left_join(
      y = profiles, by = c("Template_Gene" = "Gene")
    ) %>%
    dplyr::group_by(Sample, Serovar)
  
  logger::log_debug("Counting the capsule genes to determine repressentation of serovars.")
  kma_overview <- dplyr::summarise(
    kma_profile,
    Gene_count = sum(c(match_perfect, match_imperfect)),
    .groups = "drop_last"
  ) %>%
    # Determine which serovars are best represented relative to capsule gene counts
    dplyr::summarise(
      Serovar,
      Gene_count,
      selected = Gene_count == max(Gene_count),
      .groups = "drop"
    )
  
  logger::log_debug("Counting expected amount of capsule genes for each serovar")
  profiles_count <- dplyr::group_by(profiles, Serovar) %>%
    dplyr::summarise(capsule_count = dplyr::n())
  
  logger::log_debug("Filtering the most repressented serovar and quantifying capsule gene frequency.")
  subset(kma_overview, selected) %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(
      suggestions = dplyr::n(),
      Serovar = paste(Serovar, collapse = ","),
      count = dplyr::case_when(
        suggestions != 1 ~ NA_integer_,
        TRUE ~ unique(Gene_count)
      )
    ) %>%
    dplyr::left_join(y = profiles_count, by = "Serovar") %>%
    dplyr::summarise(
      Sample,
      Serovar = dplyr::case_when(
        is.na(count) ~ paste(Serovar, sep = ", "),
        count / capsule_count < 0.5 ~ NA_character_,
        TRUE ~ Serovar
      ),
        
      Frequency = dplyr::case_when(
        is.na(count) ~ NA_character_,
        TRUE ~ paste(count, capsule_count, sep = " of ")
      )
    )
}


resolve_genes <- function(kma_table, profiles){

  logger::log_debug("Compiling ID and Cov details for each gene.")
  gene_details <- dplyr::mutate(
    kma_table,
    gene_detailed = paste(
      Template_Gene, paste0("(", Template_Identity, "%"),
      "ID,", paste0(Template_Coverage, "%"), "COV)"
    ),
    class = stringr::str_extract(string = Template_Gene, pattern = "\\w$")
  ) %>%
    dplyr::group_by(Sample) %>%
    dplyr::arrange(class)
  
  logger::log_debug("Defining and annotating accepted gene- and partial gene-matches.")
  gene_summary <- dplyr::summarise(
    gene_details,
    Genes = dplyr::case_when(
      match_perfect ~ Template_Gene,
      match_imperfect ~ gene_detailed
    ),
    Partial_Genes = dplyr::case_when(
      match_partial ~ gene_detailed
    ),
    .groups = "keep"
  )
  
  logger::log_debug("Removing NA's and compressing accepted- and partial -genes into single columns.")
  dplyr::summarise(
    gene_summary,
    Genes = paste0(Genes[!is.na(Genes)], collapse = ", "),
    Partial_Genes = paste0(Partial_Genes[!is.na(Partial_Genes)], collapse = ", "),
    .groups = "drop"
  )

}


summarize_serovars <- function(kma_dir, serovar_config_yaml, threshold_id, threshold_cov, serovar_file){
  logger::log_info("Detecting .res files.")
  res_files <- list.files(path = kma_dir, pattern = "\\.res", full.names = TRUE, recursive = TRUE)
  
  logger::log_info("Reading and merging all .res files.")
  res_table <- purrr::map_dfr(res_files, read_res)
  kma_table <- apply_thresholds(res_table, threshold_id, threshold_cov)
  
  if (is.null(kma_table))
    stop(
      "No `.res` files detected! Check the kma results location: ",
      dirname(kma_dir) %>%
        unique
    )
  
  logger::log_info("Generating serovar profiles from profile-config file.")
  profiles <- generate_serovar_profiles(serovar_config_yaml)
  
  logger::log_info("Determining the most frequently repressented serovar")
  serovars <- resolve_serovars(kma_table, profiles)
  
  logger::log_info("Compressing genes into a presentable format")
  genes <- resolve_genes(kma_table, profiles)
  
  logger::log_info("Collecting results")
  results <- dplyr::left_join(x = serovars, y = genes, by = "Sample")
  
  if (file.exists(serovar_file)){
    logger::log_info("Reading existing results")
    results_old <- readr::read_tsv(serovar_file)
    results_merged <- dplyr::bind_rows(results, results_old)
    results <- dplyr::distinct(results_merged)
  }
  
  logger::log_info("Writing results to: ", serovar_file)
  readr::write_tsv(x = results, file = serovar_file)
  message("Success!")
}



kma_dir <- snakemake@input[["kma_dir"]]
threshold_id <- snakemake@params[["threshold_id"]]
threshold_cov <- snakemake@params[["threshold_cov"]]
serovar_file <- snakemake@output[["serovar_file"]]
dbg <- snakemake@params[["debug"]]

logger::log_threshold(level = logger::INFO)
if (dbg){
  logger::log_threshold(level = logger::DEBUG)
  
  tmp_file = file.path(
    dirname(serovar_file),
    "summarize_serovars.RData"
  )
  
  logger::log_debug("Writing snakemake R object to ", tmp_file)
  save(snakemake, file = tmp_file)
  
}

summarize_serovars(
  kma_dir = kma_dir,
  serovar_config_yaml = "config/serovar_profiles.yaml",
  threshold_id = threshold_id,
  threshold_cov = threshold_cov,
  serovar_file = serovar_file
)



