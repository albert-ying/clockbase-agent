library(pbmcapply)
library(tidyverse)
library(purrr)

Sys.setenv(VROOM_CONNECTION_SIZE = 256000)  # You can adjust the size as needed

source("src_rna/FUN.Preprocessing_for_clock_3.0.R")

mouse_rna_folder <- "./GEO_mm_rna/raw_data/"
human_rna_folder <- "./GEO_hs_rna/raw_data/"
mouse_processed_output_folder <- "./GEO_mm_rna/processed_data/"
human_processed_output_folder <- "./GEO_hs_rna/processed_data/"
# mkdir
dir.create(mouse_processed_output_folder, showWarnings = FALSE)
dir.create(human_processed_output_folder, showWarnings = FALSE)

# Load your gene list
gene_list <- read_tsv("./src_rna/Gene_list_all_4.6.txt", col_names = FALSE)$X1 |>
  as.character()

# Get the list of files in the mouse and human RNA folders
mouse_files <- list.files(mouse_rna_folder, pattern = "\\.csv\\.gz$", full.names = TRUE)
human_files <- list.files(human_rna_folder, pattern = "\\.csv\\.gz$", full.names = TRUE)

ortholog_table = read.csv("./src_rna/Table_of_orthologs.csv",header=T)
rownames(ortholog_table) = as.character(ortholog_table$Entrez.Human)

# Define the preprocessing function wrapper
gene_table_human = read_csv("./src_rna/Gene_table_human.csv") |>
  select(Ensembl, Entrez)
gene_table_mouse = read_csv("./src_rna/Gene_table_mouse.csv") |>
  select(Ensembl = ...1, Entrez)

preprocess_wrapper <- function(
  file,
  gene_list,
  output_folder,
  map_to_orthologs = FALSE,
  count_thresh = 10,
  sample_thresh = 0.2,
  map_from = "Human",
  map_to = "Mouse",
  ortholog_table,
  gene_table
) {
  gseid = gsub("\\.csv\\.gz$", "", basename(file))
  output_prefix <- file.path(output_folder, gseid)
  # Define paths for all output files
  scaled_data_path <- paste0(output_prefix, "_scaled_data.csv")
  yugene_data_path <- paste0(output_prefix, "_yugene_data.csv")
  scaled_diff_data_path <- paste0(output_prefix, "_scaled_diff_data.csv")
  yugene_diff_data_path <- paste0(output_prefix, "_yugene_diff_data.csv")
  output_folder_path <- paste0(output_folder, "/", gseid)
  correct_scaled_data_path <- paste0(output_folder_path, "/", gseid, "_scaled_data.csv")
  correct_yugene_data_path <- paste0(output_folder_path, "/", gseid, "_yugene_data.csv")
  correct_scaled_diff_data_path <- paste0(output_folder_path, "/", gseid, "_scaled_diff_data.csv")
  correct_yugene_diff_data_path <- paste0(output_folder_path, "/", gseid, "_yugene_diff_data.csv")

  # Check if all the output files already exist
  if (all(file.exists(scaled_data_path, yugene_data_path, scaled_diff_data_path, yugene_diff_data_path))) {
    message("Skipping already processed file: ", file)
    dir.create(output_folder_path, showWarnings = FALSE)
    file.rename(scaled_data_path, correct_scaled_data_path)
    file.rename(yugene_data_path, correct_yugene_data_path)
    file.rename(scaled_diff_data_path, correct_scaled_diff_data_path)
    file.rename(yugene_diff_data_path, correct_yugene_diff_data_path)
    return(NULL)
  } else if (all(file.exists(correct_scaled_data_path, correct_yugene_data_path, correct_scaled_diff_data_path, correct_yugene_diff_data_path))) {
    message("Skipping already processed file: ", file)
    return(NULL)
  } 
  dir.create(output_folder_path, showWarnings = FALSE)

  exprs_data <- read_csv(file)

  # Join expression data with gene table to get Entrez IDs
  exprs_data <- exprs_data %>%
    dplyr::left_join(gene_table, by = c("sample" = "Ensembl")) |>
    dplyr::filter(!is.na(Entrez))

  # Group by Entrez ID and summarize expression values
  exprs_data <- exprs_data %>%
    dplyr::group_by(Entrez) %>%
    dplyr::summarise(across(where(is.numeric), sum)) %>%
    tibble::column_to_rownames("Entrez")

  #Filter out low expressed genes
  r = apply(exprs_data,1,function(x)sum(x>=count_thresh)>=ncol(exprs_data)*sample_thresh)
  exprs_data = exprs_data[r,,drop=FALSE]

  #Map to mouse orthologs (for human data)
  if (map_to_orthologs==T){
    rownames(ortholog_table) = as.character(ortholog_table[,paste0("Entrez.",map_from)])
    common_genes = intersect(rownames(exprs_data),rownames(ortholog_table))
    exprs_data = exprs_data[common_genes,]
    rownames(exprs_data) = as.character(ortholog_table[common_genes,paste0("Entrez.",map_to)])
  }

  #Preprocess data for the clock
  processed_data <- preprocessing_for_clock(
    exprs_data,
    gene_list,
    log_transform = TRUE,
    show_plots = F,
    map_to_orthologs = F,
    scaling = TRUE,
    yugene = TRUE,
    differential = TRUE,
    control_group = colnames(exprs_data),
    function_path="./src_rna/"
  )

  output_prefix <- file.path(output_folder, gsub("\\.csv\\.gz$", "", basename(file)))

  # Save scaled data
  write.csv(processed_data$Scaled_data, file = correct_scaled_data_path, row.names = TRUE)

  # Save YuGene data
  write.csv(processed_data$YuGene_data, file = correct_yugene_data_path, row.names = TRUE)

  # Save scaled differential data
  write.csv(
    processed_data$Scaled_differential_data,
    file = correct_scaled_diff_data_path,
    row.names = TRUE
  )

  # Save YuGene differential data
  write.csv(
    processed_data$YuGene_differential_data,
    file = correct_yugene_diff_data_path,
    row.names = TRUE
  )
  return(NULL)
}

preprocess_wrapper = possibly(preprocess_wrapper, otherwise = ~NULL)

# Preprocess mouse RNA data
pbmclapply(
  mouse_files,
  preprocess_wrapper,
  gene_list = gene_list,
  output_folder = mouse_processed_output_folder,
  gene_table = gene_table_mouse,
  ortholog_table = ortholog_table,
  mc.cores = 10
)

# Preprocess human RNA data
pbmclapply(
  human_files,
  preprocess_wrapper,
  gene_list = gene_list,
  output_folder = human_processed_output_folder,
  map_to_orthologs = TRUE,
  gene_table = gene_table_human,
  ortholog_table = ortholog_table,
  mc.cores = 1
)
