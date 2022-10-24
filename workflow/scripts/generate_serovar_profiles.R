library(magrittr)

object_to_dataframe <- function(list_iterator, object){
  item <- object[list_iterator]
  serovar <- names(item)
  genes <- unlist(item)
  
  dplyr::tibble(Gene = genes, Serovar = serovar)
}


list_to_table <- function(object){
  purrr::map_dfr(seq_along(object), object_to_dataframe, object)
}


generate_serovar_profiles <- function(serovar_yaml = "config/serovar_profiles.yaml", profiles_file){
  serovar_profiles <- yaml::yaml.load_file(input = serovar_yaml) %>%
    list_to_table()
  
  readr::write_tsv(x = serovar_profiles, file = profiles_file)
  
  message("INFO: Done!")
}

if (snakemake@params[["debug"]]){
  file_name <- file.path(getwd(), "Serovar.RData")
  message("Storing snakemake object at: ", file_name)
  saveRDS(object = snakemake, file = file_name)
}

generate_serovar_profiles(
  profiles_file = snakemake@output[["profiles_file"]]
)
