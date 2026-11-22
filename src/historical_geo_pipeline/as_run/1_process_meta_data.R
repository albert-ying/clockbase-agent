#* metadata folders:
#* GEO_hs_dnam/metadata
#* GEO_mm_dnam/metadata
#* GEO_hs_rna/metadata
#* GEO_mm_rna/metadata

tidylog()

source('./fun/clean_metadata.R')

batch_process_file = function(files, ncore = 10) {
  pbmcapply::pbmclapply(files, function(file) {
    metadata = read_csv(file)
    processed_metadata = metadata |>
      extract_metadata()
    outfile = str_replace(file, "metadata", "processed_metadata_preharmonize")
    write_csv(processed_metadata, outfile)
    harmonized_metadata = processed_metadata |>
      harmonize_field()
    outfile = str_replace(file, "metadata", "processed_metadata")
    write_csv(harmonized_metadata, outfile)
  }, mc.cores = ncore)
}

# Test batch_process_file
# files = list.files('./GEO_hs_dnam/metadata', full.names = TRUE)
# file = paste0('./GEO_hs_dnam/metadata/', "GSE4009.csv")
# metadata = read_csv(file[[1]])
# ipaint(metadata)
# processed_metadata = metadata |>
#   extract_metadata()
# outfile = str_replace(file, "metadata", "processed_metadata")

dir.create("GEO_hs_dnam/processed_metadata", showWarnings = F)
dir.create("GEO_hs_dnam/processed_metadata_preharmonize", showWarnings = F)
metadata_files = dir("GEO_hs_dnam/metadata", full.names = T)
batch_process_file(metadata_files, ncore = 50)

dir.create("GEO_mm_dnam/processed_metadata", showWarnings = F)
dir.create("GEO_mm_dnam/processed_metadata_preharmonize", showWarnings = F)
metadata_files = dir("GEO_mm_dnam/metadata", full.names = T)
batch_process_file(metadata_files, ncore = 20)

dir.create("GEO_hs_rna/processed_metadata", showWarnings = F)
dir.create("GEO_hs_rna/processed_metadata_preharmonize", showWarnings = F)
metadata_files = dir("GEO_hs_rna/metadata", full.names = T)
batch_process_file(metadata_files, ncore = 50)

dir.create("GEO_mm_rna/processed_metadata", showWarnings = F)
dir.create("GEO_mm_rna/processed_metadata_preharmonize", showWarnings = F)
metadata_files = dir("GEO_mm_rna/metadata", full.names = T)
batch_process_file(metadata_files, ncore = 50)

## Check the number of processed metadata files ---------------------------------------------

tot_diff = c(
  setdiff(dir("GEO_hs_dnam/processed_metadata"), dir("GEO_hs_dnam/metadata")),
  setdiff(dir("GEO_mm_dnam/processed_metadata"), dir("GEO_mm_dnam/metadata")),
  setdiff(dir("GEO_hs_rna/processed_metadata"), dir("GEO_hs_rna/metadata")),
  setdiff(dir("GEO_mm_rna/processed_metadata"), dir("GEO_mm_rna/metadata"))
)

if (length(tot_diff) == 0) {
  message("All metadata files are processed!")
} else {
  message("Some metadata files are not processed!")
}
