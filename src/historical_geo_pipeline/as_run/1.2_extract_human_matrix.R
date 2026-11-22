library(geograbi)
filter = dplyr::filter

## Human DNAm --------------------

files = dir("./GEO_hs_dnam/series_matrix/", full.names = T)

out_folder = "./GEO_hs_dnam/processed_data"

if (!dir.exists(out_folder)) {
  dir.create(out_folder)
}

extract_table = function(file, out_folder) {
  gse_id = str_extract(file, "GSE[0-9]+")

  if (file.exists(glue("{out_folder}/{gse_id}.csv.gz"))) {
    print(glue("{gse_id} already exists"))
  } else {
    mat = geograbi::geograbi.read.gse.matrix(file)
    dat = mat$data |>
      as.data.frame() |>
      rownames_to_column("probe_id") |>
      as_tibble()
    write_csv(dat, glue("{out_folder}/{gse_id}.csv.gz"))
  }
}

extract_table = possibly(extract_table, otherwise = NULL)

pbmcapply::pbmclapply(
  files,
  extract_table,
  out_folder = out_folder,
  mc.cores = 2
)
