library(furrr)
plan(multicore, workers = 50)
## read in the human DNAm data -----
human_rna_age = dir("./GEO_hs_rna/clock_result", full.name = T)
mouse_rna_age = dir("./GEO_mm_rna/clock_result", full.name = T)
human_rna_meta = dir("./hs_rna_gpt-4", pattern = "dsv", full.name = T)
mouse_rna_meta = dir("./mm_rna_gpt-4", pattern = "dsv", full.name = T)

human_rna_age_tb = map_dfr(human_rna_age, read_csv, id = "file")
human_rna_meta_tb = read_delim(human_rna_meta, id = "file", col_names = F, delim = "|")

mouse_rna_age_tb = future_map_dfr(mouse_rna_age, read_csv, id = "file")
# Combine all columns named sample_id*, with the value to be the first non-NA value
mouse_rna_age_tb <- mouse_rna_age_tb %>%
  mutate(
    combined_sample_id = coalesce(!!!select(., starts_with("sample_id")))
  ) %>%
  select(-starts_with("sample_id")) |>
  select(file, sample_id = combined_sample_id, everything())

mouse_rna_meta_tb = read_delim(mouse_rna_meta, id = "file", col_names = F, delim = "|")

# write_csv(human_rna_age_tb, "result/human_rna_age.csv")
# write_csv(human_rna_meta_tb, "result/gpt4_human_rna_metadata.csv")

colnames(human_rna_meta_tb) = colnames(mouse_rna_meta_tb) = c(
  "file",
  "gsm_id",
  "PATIENT_ID",
  "race",
  "sex",
  "age",
  "genetic_info",
  "disease",
  "tissue",
  "cell_line",
  "vivo_vitro",
  "case_control",
  "group_name",
  "treatment",
  "perturbation_category"
)

human_rna_age_tb = human_rna_age_tb |>
  tidyr::extract(file, into = "gse_id", regex = "(GSE[0-9]+)")
human_rna_age_tb = human_rna_age_tb |>
  rename(gsm_id = sample_id)

human_rna_meta_tb = human_rna_meta_tb |>
  tidyr::extract(file, into = "gse_id", regex = "(GSE[0-9]+)")

mouse_rna_age_tb = mouse_rna_age_tb |>
  tidyr::extract(file, into = "gse_id", regex = "(GSE[0-9]+)")
mouse_rna_age_tb = mouse_rna_age_tb |>
  rename(gsm_id = sample_id)

mouse_rna_meta_tb = mouse_rna_meta_tb |>
  tidyr::extract(file, into = "gse_id", regex = "(GSE[0-9]+)")

write_csv(human_rna_age_tb, "result/human_rna_age.csv")
write_csv(human_rna_meta_tb, "result/gpt4_human_rna_metadata.csv")
write_csv(mouse_rna_age_tb, "result/mouse_rna_age.csv")
write_csv(mouse_rna_meta_tb, "result/gpt4_mouse_rna_metadata.csv")
