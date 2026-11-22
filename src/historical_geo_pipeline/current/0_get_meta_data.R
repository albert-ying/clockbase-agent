library(geograbi)
library(parapurrr)
tidylog()

source('./fun/clean_metadata.R')
rerun = F
## 1. Grab all datasets in GEO ---------------------------------------------------------------------
if (!file.exists('geo_fix.csv') | rerun == T) {
  datasets = geograbi.retrieve.datasets(retmax = 9000)
  #! due to the limitaion of api, the description of some datasets are truncated
  head(datasets)
  caped_des = datasets |>
    filter(str_detect(description, "more\\.\\.\\.$")) |>
    filter(organism %in% c("Homo sapiens", "Mus musculus"))
  fix_desc_ls = pbmcapply::pbmclapply(caped_des$accession, get_desc, mc.cores = 60)
  fix_desc_tb = bind_rows(fix_desc_ls)
  fix_desc_tb |>
    drop_na()
  caped_des_fix = caped_des |>
    select(-title) |>
    inner_join(fix_desc_tb, by = c("accession")) |>
    mutate(description = summary) |>
    select(-summary) |>
    drop_na(description)
    # mutate(description = get_desc(accession) |> pull(summary))
  datasets_fix = datasets |>
    filter(!accession %in% caped_des_fix$accession) |>
    bind_rows(caped_des_fix)
  write_csv(datasets_fix, "geo_fix.csv")
  datasets = datasets_fix
} else {
  datasets = read_csv("geo_fix.csv")
}

# Check if there is a overall_design column
if (! "Overall_design" %in% colnames(datasets)) {
  datasets = datasets |>
    mutate(Overall_design = NA)
}

# Check if column overall_design is empty
datasets_without_design = datasets |>
  filter(is.na(Overall_design)) |>
  pull(accession)
if (length(datasets_without_design) > 0) {
  desc_ls = pbmcapply::pbmclapply(datasets_without_design, get_desc, mc.cores = 40)
  desc_tb = bind_rows(desc_ls) |>
    select(accession, summary, Overall_design)
  datasets = datasets |>
    left_join(desc_tb, by = "accession")
  # if description is empty, use summary
  datasets = datasets |>
    mutate(description = ifelse(is.na(description), summary, description)) |>
    select(-summary)
  write_csv(datasets, "geo_fix_design.csv")
}

# paint(datasets) 2023/11/16 update
# data.frame [212592, 10]
# title       chr Single nuclei RNA-seq profiles from medulloblastoma tumors with extensive~
# description chr Medulloblastoma with extensive nodularity (MBEN) are cerebellar tumors wi~
# organism    chr Homo sapiens Homo sapiens Homo sapiens Mus musculus Gallus gallus Sacchar~
# type        chr Methylation profiling by high throughput sequencing Expression profiling ~
# platform    chr GPL21145 GPL13534 GPL33933  GPL24247 GPL3213 GPL19756 GPL24247 GPL33920 G~
# accession   chr GSE247741 GSE247736 GSE247728 GSE247727 GSE247725 GSE247583 GSE247573 GSE~
# id          chr 200247741 200247736 200247728 200247727 200247725 200247583 200247573 200~
# project     chr          
# dataset     chr          
# samples     int 9 5 NA 28 251 6 24 6 10 9

datasets |>
  separate_rows(type, sep = "; ") |>
  count(type) |>
  distinct() |> arrange(type, n) |> as.data.frame()
#                                                                type     n
# 1                                     Expression profiling by array 67442
# 2                       Expression profiling by genome tiling array   758
# 3                Expression profiling by high throughput sequencing 88625
# 4                                      Expression profiling by MPSS    24
# 5                                    Expression profiling by RT-PCR   963
# 6                                      Expression profiling by SAGE   240
# 7                                 Expression profiling by SNP array    14
# 8                       Genome binding/occupancy profiling by array   236
# 9         Genome binding/occupancy profiling by genome tiling array  2386
# 10 Genome binding/occupancy profiling by high throughput sequencing 34518
# 11                  Genome binding/occupancy profiling by SNP array    23
# 12                              Genome variation profiling by array   872
# 13                Genome variation profiling by genome tiling array  1639
# 14         Genome variation profiling by high throughput sequencing   315
# 15                          Genome variation profiling by SNP array  1531
# 16                                   Methylation profiling by array  1569
# 17                     Methylation profiling by genome tiling array  2394
# 18              Methylation profiling by high throughput sequencing  5325
# 19                               Methylation profiling by SNP array    17
# 20                                Non-coding RNA profiling by array  5535
# 21                  Non-coding RNA profiling by genome tiling array   111
# 22           Non-coding RNA profiling by high throughput sequencing  7461
# 23                                                            Other 13397
# 24                                   Protein profiling by Mass Spec     9
# 25                               Protein profiling by protein array   408
# 26                                      SNP genotyping by SNP array   889
# 27                                           Third-party reanalysis   590

## 2. Get the GEO accession numbers for each group --------------------------------
# Human DNAm datasets 3,450 rows
gse_me_hs = datasets |>
  filter(str_detect(type, "Methylation profiling by[^\\;]*array") & organism == "Homo sapiens")
write_csv(gse_me_hs, "gse_me_hs.csv")

# Mouse DNAm datasets 1,960 rows
gse_me_mm = datasets |>
  filter(str_detect(type, "Methylation profiling by high throughput sequencing") & organism == "Mus musculus")
write_csv(gse_me_mm, "gse_me_mm.csv")

# Human RNA-seq datasets 59,503 rows
gse_rna_hs = datasets |>
  filter(str_detect(type, "Expression profiling by [array|high throughput sequencing]") & organism == "Homo sapiens")
write_csv(gse_rna_hs, "gse_rna_hs.csv")

# mouse RNA-seq datasets 51,771 rows
gse_rna_mm = datasets |>
  filter(str_detect(type, "Expression profiling by [array|high throughput sequencing]") & organism == "Mus musculus")
write_csv(gse_rna_mm, "gse_rna_mm.csv")

## 3. Download the series file ------------------------------------------------------------------
filenames = geograbi.download.series.files(path = "./GEO_hs_dnam/series_matrix", gses = gse_me_hs$accession)
filenames = geograbi.download.series.files(path = "./GEO_mm_dnam/series_matrix", gses = gse_me_mm$accession)
filenames = geograbi.download.series.files(path = "./GEO_hs_rnaseq/series_matrix", gses = gse_rna_hs$accession)
filenames = geograbi.download.series.files(path = "./GEO_mm_rnaseq/series_matrix", gses = gse_rna_mm$accession)

## Parse the metadata ---------------------------------------------------------------------------
# Define the metadata parser
save_meta = function(matrix_file, dest) {
  geo = basename(matrix_file) |> str_extract("GSE[0-9]+")
  dest_file = file.path(dest, paste0(geo, ".csv"))
  if (file.exists(dest_file)) {
    message("File ", dest_file, " exists, skip")
    return()
  } else {
    message("Parsing ", matrix_file)
    meta = geograbi.get.samples(matrix_file)
    meta_tb = meta |>
      rownames_to_column("GSM") |>
      as_tibble()
    message("Writing to ", dest_file)
    write_csv(meta_tb, dest_file)
  }
}

# Save the metadata
geo_hs_dnam = dir("GEO_hs_dnam/series_matrix", full.names = T)
walk(geo_hs_dnam, save_meta, dest = "GEO_hs_dnam/metadata")

geo_mm_dnam = dir("GEO_mm_dnam/series_matrix", full.names = T)
walk(geo_mm_dnam, save_meta, dest = "GEO_mm_dnam/metadata")

geo_hs_rna = dir("GEO_hs_rna/series_matrix", full.names = T)
walk(geo_hs_rna, save_meta, dest = "GEO_hs_rna/metadata")

geo_mm_rna = dir("GEO_mm_rna/series_matrix", full.names = T)
walk(geo_mm_rna, save_meta, dest = "GEO_mm_rna/metadata")


# # cleanup
# valid_gse = gse_rna_hs$accession
# for (file in dir("./GEO_hs_rna/series_matrix/", full.name = T)) {
#   gse = file |> str_extract("GSE[0-9]*")
#   if (! gse %in% valid_gse) {
#     message(gse)
#     file.remove(file)
#   }
# }

# valid_gse = gse_rna_mm$accession
# for (file in dir("./GEO_mm_rna/series_matrix/", full.name = T)) {
#   gse = file |> str_extract("GSE[0-9]*")
#   if (! gse %in% valid_gse) {
#     message(gse)
#     file.remove(file)
#   }
# }

# length(dir("./GEO_mm_rna/series_matrix/", full.name = T))

# length(dir("./GEO_hs_rna/series_matrix/", full.name = T))
