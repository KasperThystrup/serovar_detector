library(magrittr)


read_res <- function(res_file, dbg){
  
  logger::log_debug("Performing pattern matching to determine sample name from file: ", res_file)
  file_matches <- res_file %>%
    basename %>%
    stringr::str_match_all(pattern = "^(?<sample>[A-Za-z0-9_-]+)\\.(?<ext>\\w+)$") %>%
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
  
  dplyr::mutate(
    res_table,
    match_perfect = Template_Identity == 100 & Template_Coverage == 100,
    match_imperfect = !match_perfect & (Template_Identity >= threshold_id & Template_Coverage >= threshold_cov), # Manual thresholds!!
    match_partial = !match_imperfect & !match_perfect
  )
}


summarise_serovars <- function(kma_table, profiles, dbg){
  
  logger::log_debug("Merging kma results with serovar profiles table.")
  kma_profile <- dplyr::select(
    kma_table, 
    Sample, Template_Gene, Template_Strain, Template_Serovar
  ) %>%
    dplyr::left_join(
      x = kma_table, y = profiles, by = c("Template_Gene" = "Gene")
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
      serovar_count = max(Gene_count),
      selected = Gene_count == serovar_count,
      .groups = "drop"
    )
  
  logger::log_debug("Counting expected amount of capsule genes for each serovar")
  profiles_count <- dplyr::group_by(profiles, Serovar) %>%
    dplyr::summarise(capsule_count = dplyr::n())
  if (dbg)
    message(head(profiles_count))
  
  logger::log_debug("Filtering the most repressented serovar and quantifying capsule gene frequency.")
  subset(kma_overview, selected) %>%
    dplyr::select(Sample, Serovar, serovar_count) %>%
    dplyr::left_join(y = profiles_count, by = "Serovar") %>%
    dplyr::summarise(
      Sample,
      Serovar = paste(Serovar, sep = ", "),
      Frequency = paste(serovar_count, capsule_count, sep = " of ")
    )
}


summarise_genes <- function(kma_table, profiles, dbg){

  logger::log_debug("Compiling ID and Cov details for each gene.")
  gene_details <- dplyr::mutate(
    kma_table,
    gene_detailed = paste(
      Template_Gene, paste0("(", Template_Identity, "%"),
      "ID,", paste0(Template_Coverage, "%"), "COV)"
    )
  ) %>%
    dplyr::group_by(Sample)
  
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


summarize_serovars <- function(kma_dir, profiles_file, serovar_file, threshold_id, threshold_cov, dbg){
  if (dbg){
    logger::log_threshold(logger::DEBUG)
    file_name <- file.path(getwd(), "summarize.RData")
    logger::log_debug("Storing snakemake object at: ", file_name)
    save(snakemake, file = file_name)
  }
  
  logger::log_info("Detecting .res files.")
  res_files <- list.files(path = kma_dir, pattern = "\\.res", full.names = TRUE, recursive = TRUE)
  
  logger::log_info("Reading and merging all .res files.")
  res_table <- purrr::map_dfr(res_files, read_res, dbg)
  kma_table <- apply_thresholds(res_table, threshold_id, threshold_cov)
  
  logger::log_info("Importing serovar profiles from file: ", profiles_file)
  profiles <- readr::read_tsv(file = profiles_file, col_types = readr::cols())
  
  logger::log_info("Determining the most frequently repressented serovar")
  serovars <- summarise_serovars(kma_table, profiles, dbg)
  
  logger::log_info("Compressing genes into a presentable format")
  genes <- summarise_genes(kma_table, profiles, dbg)
  
  logger::log_info("Collecting results")
  results <- dplyr::left_join(x = serovars, y = genes, by = "Sample")
  message("Success!")
  
  if (file.exists(serovar_file)){
    logger::log_info("Reading existing results")
    results_old <- readr::read_tsv(serovar_file)
    results_merged <- dplyr::bind_rows(results, results_old)
    results <- dplyr::distinct(results_merged)
  }
  
  logger::log_info("Writing results to: ", serovar_file)
  readr::write_tsv(x = results, file = serovar_file)
  message("Done!")
}


summarize_serovars(
  kma_dir = snakemake@input[["kma_dir"]],
  profiles_file = snakemake@input[["profiles_file"]],
  threshold_id = snakemake@params[["threshold_id"]],
  threshold_cov = snakemake@params[["threshold_cov"]],
  serovar_file = snakemake@output[["serovar_file"]],
  dbg = snakemake@params[["debug"]]
)



