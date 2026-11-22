library(pbmcapply)
tidylog()
source("./fun/supplementary_data_parsing.R")
source("./fun_rna/MASTER.RNAseq_preprocessing.R")

## RNAs ----------------

## Function for supplementary file processing 

RNA_process_wrapper = function(supp_folder, out_folder, gene_ref = gene_ref_mm) {
  supp_files = dir(supp_folder, full.names = T)
  gseid = str_match(supp_folder, "GSE[0-9]*") |> pluck(1)
  gsmid = str_match(supp_files, "GSM[0-9]*") |> as.vector()

  parsed_expr = map2_dfr(supp_files, gsmid, ~parse_rna_count(.x, gene_ref) |>
    mutate(gsmid = .y))

  count_table = pivot_wider(
    parsed_expr,
    names_from = gsmid,
    values_from = count,
    values_fn = sum
  ) |>
    as.data.frame() |>
    column_to_rownames("gene")

  normed_expr = RNAseq_exprs_prep(count_table)
  # Need to remove NAs
  scaled = t(normed_expr[[1]]) |>
    as_tibble() |>
    mutate(gene = colnames(normed_expr[[1]]), .before = 1) |>
    arrange(gene) |>
    drop_na(gene) |>
    group_by(gene) |>
    mutate_at(vars(-gene), sum) |>
    distinct()
  normal = t(normed_expr[[2]]) |>
    as_tibble() |>
    mutate(gene = colnames(normed_expr[[2]]), .before = 1) |>
    arrange(gene) |>
    drop_na(gene) |>
    group_by(gene) |>
    mutate_at(vars(-gene), sum) |>
    distinct()
  write_csv(scaled, file.path(out_folder, glue("{gseid}_scaled.csv")))
  write_csv(normal, file.path(out_folder, glue("{gseid}_normal.csv")))
}

possible_RNA_process_wrapper = possibly(RNA_process_wrapper, otherwise = NULL)

safe_RNA_process_wrapper = safely(RNA_process_wrapper)
Sys.setenv(VROOM_CONNECTION_SIZE = 2^20)

## Mouse RNA-seq datasets --------
folder = "GEO_mm_rna"
gene_ref_mm = read_csv("./fun_rna/Gene_table_mouse.csv")
untar_folder = glue("{folder}/supplementary_file_untar")
supp_folders = dir(untar_folder, full.names = T)
out_folder = "GEO_mm_rna/processed_data"

done_file = dir(out_folder) |>
  str_extract("GSE[0-9]+") |>
  unique()

length(dir(untar_folder))
length(done_file)

# supp_folder = supp_folders[which(str_detect(supp_folders, "GSE100237"))]
for (supp_folder in supp_folders) {
  gse_id = str_extract(supp_folder, "GSE[0-9]+")
  print(gse_id)
  if (gse_id %in% done_file) {
    next
  } else {
    error = safe_RNA_process_wrapper(supp_folder, out_folder, gene_ref = gene_ref_mm)
    error = error[[2]] |> as.character()
    fileConn <- file(glue("{out_folder}/{str_extract(supp_folder, 'GSE[0-9]+')}.log"))
    writeLines("LOG", fileConn)
    writeLines(error, fileConn)
    close(fileConn)
  }
}
# pbmclapply(supp_folders, function(supp_folder) {
#   gse_id = str_extract(supp_folder, "GSE[0-9]+")
#   print(gse_id)
#   if (gse_id %in% done_file) {
#     return(NULL)
#   } else {
#     error = safe_RNA_process_wrapper(supp_folder, out_folder, gene_ref = gene_ref_mm)
#     error = error[[2]] |> as.character()
#     fileConn <- file(glue("{out_folder}/{str_extract(supp_folder, 'GSE[0-9]+')}.log"))
#     writeLines("LOG", fileConn)
#     writeLines(error, fileConn)
#     close(fileConn)
#   }
#   return(NULL)
# }, mc.cores = 10)

## Human RNA-seq datasets -----------
folder = "GEO_hs_rna"
gene_ref_hs = read_csv("./fun_rna/Gene_table_human.csv")
untar_folder = glue("{folder}/supplementary_file_untar")
supp_folders = dir(untar_folder, full.names = T)
out_folder = "GEO_hs_rna/processed_data"

done_file = dir(out_folder) |>
  str_extract("GSE[0-9]+") |>
  unique()
# all_ext = supp_folders |>
#   map(~dir(.x)[1]) |>
#   unlist() |>
#   str_remove(".gz") |>
#   str_extract("\\.[a-z]*$")
# tibble(all_ext) |>
#   count(all_ext) |>
#   arrange(desc(n))
supp_folder = supp_folders[which(str_detect(supp_folders,"GSE139597"))]
# supp_folder = supp_folders[101]
for (supp_folder in supp_folders) {
  gse_id = str_extract(supp_folder, "GSE[0-9]+")
  print(gse_id)
  if (gse_id %in% done_file) {
    next
  } else {
    # debug safe_RNA_process_wrapper

    error = safe_RNA_process_wrapper(supp_folder, out_folder, gene_ref = gene_ref_hs)
    error = error[[2]] |> as.character()
    fileConn = file(glue("{out_folder}/{str_extract(supp_folder, 'GSE[0-9]+')}.log"))
    writeLines("LOG", fileConn)
    writeLines(error, fileConn)
    close(fileConn)
  }
}
# pbmclapply(supp_folders, function(supp_folder) {
#   gse_id = str_extract(supp_folder, "GSE[0-9]+")
#   print(gse_id)
#   if (gse_id %in% done_file) {
#     next
#   } else {
#     error = safe_RNA_process_wrapper(supp_folder, out_folder, gene_ref = gene_ref_hs)
#     error = error[[2]] |> as.character()
#     fileConn <- file(glue("{out_folder}/{str_extract(supp_folder, 'GSE[0-9]+')}.log"))
#     writeLines("LOG", fileConn)
#     writeLines(error, fileConn)
#     close(fileConn)
#   }
# }, mc.cores = 10)
## Mouse DNAm datasets -------------
