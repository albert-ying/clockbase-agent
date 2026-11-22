library(fuzzyjoin)
library(ClockBasis)
library(pbmcapply)
tidylog()

## Load GEO data
geo_folder = "../02_GEO_DNAm/GEO_mouse_processed"
geo_ids = dir(geo_folder) |>
  str_extract("GSE\\d+") |>
  unique()

## Calculate age for each sample
out_dir = "./temp/GEO_mm_dnam/feature_table/"
dir.create(out_dir, recursive = T)
# for (i in 1:length(geo_ids)) {
save_mm_feature_tb = function(i) {
  # ? test setting start ----
  # i = 1
  # ? test setting end ----
  id = geo_ids[i]
  print(id)
  gsm_files = dir(file.path(geo_folder, id), full.names = T)

  sample_ls = map(gsm_files, ~{
    gsm_id = str_extract(basename(.x), "GSM\\d+")
    read_csv(.x, col_names = c("probe_id", gsm_id))
  })
  sample_tb = read_csv(gsm_files, id = "gsm_id")
  sample_tb_processed = sample_tb |>
    mutate(gsm_id = str_extract(basename(gsm_id), "GSM\\d*")) |>
    pivot_wider(names_from = gsm_id, values_from = level) |>
    select(chr_bp, everything())
  write_csv(sample_tb_processed, paste0(out_dir, id, ".csv"))
  return(NULL)
}
save_mm_feature_tb = possibly(save_mm_feature_tb, otherwise = NULL)
for (i in seq_along(geo_ids)) {
  save_mm_feature_tb(i)
  message(i / length(geo_ids))
}
pbmclapply(seq_along(geo_ids), save_mm_feature_tb, mc.cores = 5)
