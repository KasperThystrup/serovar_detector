load(srv_2_z("/srv/data/AS/ResGroup/Tools/serovar_detector/mini/summarize_serovars.RData"))

assembly_files <- snakemake@input[["assembly_results"]] # snakemake@input[["assembly_files"]]
reads_files <- snakemake@input[["reads_results"]] # snakemake@input[["reads_files"]]
threshold <- snakemake@params[["threshold"]]
debug <- snakemake@params[["debug"]]
results_wide_file <- snakemake@output[["results_wide_file"]]
results_long_file <- snakemake@output[["results_long_file"]]


import_results <- function(results_file){
  # Read results file and repair column names
  results_raw <- readr::read_tsv(
    file = results_file,
    col_types = "ciiiddddddd",
    name_repair = ~ janitor::make_clean_names(string = .x, replace = c("\\#" = ""))
  )
  
  # Determine sample id
  sample_id <- basename(results_file) |>
    tools::file_path_sans_ext()
  
  # Prepare imported data
  results <- tidyr::separate(
    data = dplyr::mutate(results_raw, sample_id = sample_id, .before = 1),
    col = template, into = c("gene", "strain", "feature"), sep = "_"
  )
  
  # Remove feature column
  dplyr::select(results, -feature) ## feature is a helper string for the database
}


filter_results <- function(results){
  # Filter for top score by 1: Coverage, 2: Identity, and 3: Score
  results_top <- dplyr::filter(
    .data = dplyr::group_by(results, sample_id, gene),
    template_coverage == max(template_coverage)
  ) |>
    dplyr::filter(template_identity == max(template_identity)) |>
    dplyr::filter(score == max(score))
  
  # Arrange top results, filter by threshold and select top hit only
  dplyr::filter(
    dplyr::arrange(results_top, template_coverage, template_identity, score),
    template_coverage >= threshold & template_coverage >= threshold
  ) |>
    dplyr::slice_head(n = 1)
}


implement_profiles <- function(profile_file, profile){
  # Read profiles file directly into a data frame
  profiles_raw <- yaml::read_yaml(profile_file) |>
    plyr::ldply(tibble::as_tibble)
  
  # Generate clean column names
  profiles_clean <- dplyr::select(
    .data = dplyr::tibble(profiles_raw),
    "gene" = value,
    "feature" = .id
  )
  
  # Add feature to the frame
  dplyr::mutate(profiles_clean, profile = profile)
}


merge_profiles <- function(results, serovar_profiles){
  # Merge results and profiles
  results_serovar <- dplyr::left_join(
    x = results, y = serovar_profiles, by = "gene",
    relationship = "many-to-many"
  )
  
  results_merged <- results_serovar # dplyr::left_join next profile
  dplyr::relocate(results_merged, feature, profile, .before = score) 
}

results_assemblies <- purrr::map_dfr(.x = srv_2_z(assembly_files), .f = import_results)
results_reads <- purrr::map_dfr(.x = srv_2_z(reads_files), .f = import_results)
results <- dplyr::bind_rows(results_assemblies, results_reads)

results_filtered <- filter_results(results)

serovar_profiles <- implement_profiles(profile_file = profile_file, profile = "Serovar")

results_long <- merge_profiles(results_filtered, serovar_profiles)


## Long table done -> Export to results_long_file
### Make gene label (%ID and %Cov) after merging? (allready grouped per gene)
### Make summary table now :-D

profiles_members <- dplyr::summarise(
  .data = dplyr::group_by(serovar_profiles, profile, feature),
  members = dplyr::n(),
  .groups = "keep"
)

results_count <- dplyr::summarise(
  .data = dplyr::group_by(results_long, sample_id, profile, feature),
  count = dplyr::n(),
  .groups = "keep"
)

results_top <- dplyr::mutate(
  .data = dplyr::left_join(results_count, profiles_members),
  freq = count/members
)

dplyr::filter(results_top, freq == max(freq))

results_summary <- dplyr::summarise(
  .data = dplyr::group_by(results_members, sample_id, profile, feature),
  frequency = glue::glue("{dplyr::n()}/{members}")
)
