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
  output_prefix <- file.path(output_folder, gsub("\\.csv\\.gz$", "", basename(file)))
  # Define paths for all output files
  scaled_data_path <- paste0(output_prefix, "_scaled_data.csv")
  yugene_data_path <- paste0(output_prefix, "_yugene_data.csv")
  scaled_diff_data_path <- paste0(output_prefix, "_scaled_diff_data.csv")
  yugene_diff_data_path <- paste0(output_prefix, "_yugene_diff_data.csv")

  # Check if all the output files already exist
  if (all(file.exists(scaled_data_path, yugene_data_path, scaled_diff_data_path, yugene_diff_data_path))) {
    message("Skipping already processed file: ", file)
    return(NULL)
  }

  exprs_data <- read_csv(file)

  # Join expression data with gene table to get Entrez IDs
  exprs_data <- exprs_data %>%
    dplyr::left_join(gene_table, by = c("gene" = "Ensembl")) |>
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
  processed_data <- preprocessing_for_clock(exprs_data, gene_list, log_transform = TRUE, show_plots = F, 
                                            map_to_orthologs = F,
                                            scaling = TRUE, yugene = TRUE, differential = TRUE,
                                            control_group = colnames(exprs_data),
                                            function_path="./src_rna/")

  output_prefix <- file.path(output_folder, gsub("\\.csv\\.gz$", "", basename(file)))

  # Save scaled data
  write.csv(processed_data$Scaled_data, file = paste0(output_prefix, "_scaled_data.csv"), row.names = TRUE)

  # Save YuGene data
  write.csv(processed_data$YuGene_data, file = paste0(output_prefix, "_yugene_data.csv"), row.names = TRUE)

  # Save scaled differential data
  write.csv(
    processed_data$Scaled_differential_data,
    file = paste0(output_prefix, "_scaled_diff_data.csv"),
    row.names = TRUE
  )

  # Save YuGene differential data
  write.csv(
    processed_data$YuGene_differential_data,
    file = paste0(output_prefix, "_yugene_diff_data.csv"),
    row.names = TRUE
  )
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
  mc.cores = 5
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
  mc.cores = 5
)

# find the files for human that failed in preprocess

mouse_gse = str_extract(mouse_files, "GSE[0-9]+") |> unique()
processed_mouse_gse = mouse_processed_output_folder |> dir() |> str_extract("GSE[0-9]+") |> unique()
failed_mouse_gse = setdiff(mouse_gse, processed_mouse_gse)

example1 = glue::glue("./GEO_mm_rna/raw_data/{failed_mouse_gse[10000]}.csv.gz")
example2 = glue::glue("./GEO_mm_rna/raw_data/{processed_mouse_gse[100]}.csv.gz")

preprocess_wrapper(
  file = example2,
  gene_list = gene_list,
  output_folder = mouse_processed_output_folder,
  map_to_orthologs = FALSE,
  gene_table = gene_table_mouse,
  ortholog_table = ortholog_table
)

human_gse = str_extract(human_files, "GSE[0-9]+") |> unique()
processed_gse = human_processed_output_folder |> dir() |> str_extract("GSE[0-9]+") |> unique()

failed_gse = setdiff(human_gse, processed_gse)

example1 = glue::glue("./GEO_hs_rna/raw_data/{failed_gse[1000]}.csv.gz")
example2 = glue::glue("./GEO_hs_rna/raw_data/{processed_gse[1]}.csv.gz")



# run one human example to debug
preprocess_wrapper(
  file = example1,
  gene_list = gene_list,
  output_folder = human_processed_output_folder,
  map_to_orthologs = TRUE,
  gene_table = gene_table_human,
  ortholog_table = ortholog_table
)

##

# Automatically find file paths within output folders for checking the preprocessing results using dir command
mouse_scaled_data_files <- dir(path = mouse_processed_output_folder, pattern = "scaled_data.csv", full.names = TRUE)
mouse_yugene_data_files <- dir(path = mouse_processed_output_folder, pattern = "yugene_data.csv", full.names = TRUE)
mouse_scaled_diff_data_files <- dir(path = mouse_processed_output_folder, pattern = "scaled_diff_data.csv", full.names = TRUE)
mouse_yugene_diff_data_files <- dir(path = mouse_processed_output_folder, pattern = "yugene_diff_data.csv", full.names = TRUE)

human_scaled_data_files <- dir(path = human_processed_output_folder, pattern = "scaled_data.csv", full.names = TRUE)
human_yugene_data_files <- dir(path = human_processed_output_folder, pattern = "yugene_data.csv", full.names = TRUE)
human_scaled_diff_data_files <- dir(path = human_processed_output_folder, pattern = "scaled_diff_data.csv", full.names = TRUE)
human_yugene_diff_data_files <- dir(path = human_processed_output_folder, pattern = "yugene_diff_data.csv", full.names = TRUE)

# Assuming we want to check the first file from each output type for mouse data
check_mouse_scaled_data <- read.csv(mouse_scaled_data_files[1], row.names = 1)
print(head(check_mouse_scaled_data))

check_mouse_yugene_data <- read.csv(mouse_yugene_data_files[1], row.names = 1)
print(head(check_mouse_yugene_data))

check_mouse_scaled_diff_data <- read.csv(mouse_scaled_diff_data_files[1], row.names = 1)
print(head(check_mouse_scaled_diff_data))

check_mouse_yugene_diff_data <- read.csv(mouse_yugene_diff_data_files[1], row.names = 1)
print(head(check_mouse_yugene_diff_data))

# Assuming we want to check the first file from each output type for human data
check_human_scaled_data <- read.csv(human_scaled_data_files[1], row.names = 1)
print(head(check_human_scaled_data))

check_human_yugene_data <- read.csv(human_yugene_data_files[1], row.names = 1)
print(head(check_human_yugene_data))

check_human_scaled_diff_data <- read.csv(human_scaled_diff_data_files[1], row.names = 1)
print(head(check_human_scaled_diff_data))

check_human_yugene_diff_data <- read.csv(human_yugene_diff_data_files[1], row.names = 1)
print(head(check_human_yugene_diff_data))
