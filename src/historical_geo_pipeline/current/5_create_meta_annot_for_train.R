#* metadata folders:
#* GEO_hs_dnam/metadata
#* GEO_mm_dnam/metadata
#* GEO_hs_rna/metadata
#* GEO_mm_rna/metadata

tidylog()

source('./fun/clean_metadata.R')


extract_metadata()
batch_process_file = function(files, ncore = 10) {
  pbmcapply::pbmclapply(files, function(file) {
    metadata = read_csv(file)
    processed_metadata = metadata |>
      extract_metadata() |>
      harmonize_field()
    outfile = str_replace(file, "metadata", "processed_metadata")
    write_csv(processed_metadata, outfile)
  }, mc.cores = ncore)
}

dir.create("GEO_hs_dnam/metadata_training", showWarnings = F)

all_meta_file = dir("GEO_hs_dnam/metadata", full.names = T)
all_processed_meta_file = dir("GEO_hs_dnam/processed_metadata", full.names = T)
all_geo = str_extract(all_meta_file, "GSE[0-9]+")

geo
which(all_geo == geo)
for (i in 1374:length(all_geo)) {
  geo = all_geo[i]
  meta = read_csv(str_c("GEO_hs_dnam/metadata/", geo, ".csv")) |>
    extract_metadata()
  processed_meta = read_csv(str_c("GEO_hs_dnam/processed_metadata/", geo, ".csv"))
  target_meta = processed_meta |>
    mutate(
      GSM,
      sex = harmonized.sex,
      race = harmonized.race,
      age = harmonized.age,
      disease = harmonized.disease,
      genetic_info = harmonized.genetic_information,
      cell_line = harmonized.cell_line,
      tissue = source_name,
      vivo_vitro = NA,
      treatment = harmonized.treatment,
      case_control = NA,
      perturbation_category = NA
    ) |>
    select(GSM, sex, race, age, disease, genetic_info, cell_line,
      tissue, vivo_vitro, treatment, case_control, perturbation_category)
  desc = get_desc(geo)
  write_csv(desc, str_c("GEO_hs_dnam/metadata_training/", geo, "_desc.csv"))
  write_csv(meta, str_c("GEO_hs_dnam/metadata_training/", geo, "_meta.csv"))
  write_csv(target_meta, str_c("GEO_hs_dnam/metadata_training/", geo, "_target.csv"))
}

