library(magrittr)


read_res <- function(res_file){
  
  logger::log_debug("Performing pattern matching to determine sample name from file: ", res_file)
  file_matches <- res_file %>%
    basename %>%
    stringr::str_match_all(pattern = "^(?<sample>[^\\.]+(\\.[^\\.]+)*)\\.(?<ext>\\w+)$") %>%
    as.data.frame()

  sample_name <- file_matches$sample
  
  logger::log_debug("Reading and updating headers.")
  header <- res_file %>%
    readr::read_lines(n_max = 1) %>%
    strsplit(split = "\t") %>%
    unlist
  
  header[1] <- "Template"
  
  logger::log_debug("Determining file modification date.")
  timestamp <- file.mtime(res_file)
  mod_date <- format(timestamp, format = "%Y-%m-%d")
  
  logger::log_debug("Reading file.\n\n")
  readr::read_tsv(file = res_file, col_names = header, col_types = "ciiiddddddd", skip = 1) %>%
    tidyr::separate(
      col = Template,
      into = c("Template_Gene", "Template_Strain", "Template_Serovar"),
      sep = "_"
    ) %>%
    # Adding sample name column
    dplyr::mutate(Sample = sample_name, Date = mod_date)
}


apply_thresholds <- function(res_table, threshold){
  logger::log_debug("Applying thresholds")
  
  empty_table <- nrow(res_table) == 0
  
  if (empty_table)
    return(NULL)
  
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


resolve_serovars <- function(kma_table, profiles){
  
  logger::log_debug("Extracting Date information")
  kma_dates <- dplyr::select(kma_table, Sample, Date) %>%
    dplyr::distinct()
  
  logger::log_debug("Merging kma results with serovar profiles table.") ### Warnings!!!
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
    dplyr::reframe(
      Serovar,
      Gene_count,
      selected = Gene_count == max(Gene_count),
      .groups = "drop"
    )
  
  logger::log_debug("Counting expected amount of capsule genes for each serovar")
  profiles_count <- dplyr::group_by(profiles, Serovar) %>%
    dplyr::summarise(capsule_count = dplyr::n(), .groups = "keep") ## .groups added
  
  logger::log_debug("Filtering the most repressented serovar and quantifying capsule gene frequency.")
  serovar_suggestions <- subset(kma_overview, selected) %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(
      suggestions = dplyr::n(),
      Serovar = paste(Serovar, collapse = ","),
      count = dplyr::case_when(
        suggestions != 1 ~ NA_integer_,
        TRUE ~ unique(Gene_count)
      ),
      .groups = "keep"
    )
  
  serovars_raw <- dplyr::left_join(
    x = serovar_suggestions, y = profiles_count, by = "Serovar"
  ) %>%
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
  
  serovars <- dplyr::inner_join(x = serovars_raw, kma_dates, by = "Sample") %>%
    dplyr::relocate(Date, .after = Sample)
  
  kma_merged <- dplyr::left_join(kma_profile, serovars, by = "Sample") %>%
    dplyr::group_by(Sample, Template_Gene)
  
  kma_detailed <- dplyr::mutate(
    kma_merged,
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
  ) %>%
    dplyr::group_by(Sample) %>%
    dplyr::arrange(class)
  
  logger::log_debug("Defining and annotating accepted gene- and partial gene-matches.")
  kma_genes <- dplyr::reframe(
    kma_detailed,
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
  ) %>%
    dplyr::group_by(Sample)
  
  logger::log_debug("Removing NA's and compressing accepted- and partial -genes into single columns.")
  genes <- dplyr::summarise(
    kma_genes,
    Serovar_match = paste0(unique(Serovar_match[!is.na(Serovar_match)]), collapse = ", "),
    Serovar_partial = paste0(unique(Serovar_partial[!is.na(Serovar_partial)]), collapse = ", "),
    Others_match = paste0(unique(Others_match[!is.na(Others_match)]), collapse = ", "),
    Others_partial = paste0(unique(Others_partial[!is.na(Others_partial)]), collapse = ", "),
    .groups = "drop"
  )
  
  dplyr::left_join(serovars, genes, by = "Sample")
    
}
  

summarize_serovars <- function(kma_files, serovar_config_yaml, threshold, blacklisting, blacklist_file, serovar_file){
  logger::log_info("Reading and merging all .res files.")
  res_table <- purrr::map_dfr(kma_files, read_res)
  kma_table <- apply_thresholds(res_table, threshold)
  
  if (is.null(kma_table))
    stop(
      "No `.res` files detected! Check the kma results location: ",
      dirname(kma_dir) %>%
        unique
    )
  
  logger::log_info("Generating serovar profiles from profile-config file.")
  profiles <- generate_serovar_profiles(serovar_config_yaml)
  
  logger::log_info("Determining the most frequently repressented serovars and serovar genes.")
  results <- resolve_serovars(kma_table, profiles)
  
  if (blacklisting){
    blacklisted_samples <- dplyr::select(results, Sample)
    
    if (file.exists(blacklist_file)){
      blacklist_input <- readr::read_tsv(blacklist_file)
    
      blakclisted_samples <- dplyr::bind_rows(blacklist_input, blacklisted_samples)
      
      any_blacklist_duplicates <- dplyr::pull(blacklisted_samples, var = Sample) %>%
        duplicated %>%
        any
      
      if (any_blacklist_duplicates){
        warning("Duplicated blacklist sampels detected. Please check the `config/blacklist.txt`")
        blacklisted_samples <- dplyr::distinct(blacklisted_samples)
    }
    
    readr::write_tsv(x = blacklisted_samples, file = blacklist_file)
  }
  if (file.exists(serovar_file)){
    logger::log_info("Reading existing results")
    results_old <- readr::read_tsv(serovar_file)
    old_samples <- dplyr::pull(results_old, Sample)
    results_new <- subset(results, !(Sample %in% old_samples))
    results_merged <- dplyr::bind_rows(results_new, results_old)
  }
  
  logger::log_info("Writing results to: ", serovar_file)
  readr::write_tsv(x = dplyr::arrange(results, Date, Sample), file = serovar_file)
  message("Success!")
}

assembly_results <- snakemake@input[["assembly_results"]]
reads_results <- snakemake@input[["reads_results"]]
threshold <- snakemake@params[["threshold"]]
blacklisting <- snakemake@params[["blacklisting"]]
serovar_file <- snakemake@output[["serovar_file"]]
dbg <- snakemake@params[["debug"]]

kma_files <- c(assembly_results, reads_results)

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
  kma_files = kma_files,
  serovar_config_yaml = "config/serovar_profiles.yaml",
  threshold = threshold,
  blacklisting = blacklisting,
  blacklist_file = "config/blacklist.tsv",
  serovar_file = serovar_file
)
