library(ClockBasis)
library(methylCIPHER)
library(geograbi)
filter = dplyr::filter

## Human DNAm --------------------

files = dir("./GEO_hs_dnam/series_matrix/", full.names = T)

out_folder = "./GEO_hs_dnam/sample_age"

if (!dir.exists(out_folder)) {
  dir.create(out_folder)
}

run_human_clock = function(file, out_folder) {
  gse_id = str_extract(file, "GSE[0-9]+")

  if (file.exists(glue("{out_folder}/{gse_id}.csv"))) {
    print(glue("{gse_id} already exists"))
  } else {
    mat = geograbi::geograbi.read.gse.matrix(file)
    if (is.null(dim(mat$data)) || nrow(mat$data) == 0L || ncol(mat$data) == 0L) {
      stop(glue("AUDIT_FAIL: {gse_id} has no measured features or samples; clock execution refused"))
    }
    if (any(colSums(is.finite(as.matrix(mat$data))) == 0L)) {
      stop(glue("AUDIT_FAIL: {gse_id} contains samples without observed features"))
    }

    dat = mat$data |>
      as.data.frame() |>
      rownames_to_column("probe_id") |>
      as_tibble() |>
      ClockBasis::do_dnam_clock_human(ncores = 5)

    write_csv(dat, glue("{out_folder}/{gse_id}.csv"))
  }
}

possible_run_human_clock = possibly(run_human_clock, otherwise = NULL, quiet = FALSE)

pbmcapply::pbmclapply(
  files,
  possible_run_human_clock,
  out_folder = out_folder,
  mc.cores = 2
)

## Mouse DNAm --------------------

# files = dir("./GEO_mm_dnam/processed_data/", full.names = T)

# out_folder = "./GEO_mm_dnam/sample_age"

# if (!dir.exists(out_folder)) {
#   dir.create(out_folder)
# }

# file = files[1]
# run_mouse_clock = function(file, out_folder) {
#   gse_id = str_extract(file, "GSE[0-9]+")
#   mat = read_csv(file)
#   dat = mat |>
#     ClockBasis::do_dnam_clock_mouse()
#   write_csv(dat, glue("{out_folder}/{gse_id}.csv"))
# }

# possible_run_mouse_clock = possibly(run_mouse_clock, otherwise = NULL)

# # for (file in files) {
# #   which(file == files) |>
# #     print()
# #   possible_run_human_clock(file, out_folder)
# # }

# pbmcapply::pbmclapply(
#   files,
#   possible_run_mouse_clock,
#   out_folder = out_folder,
#   mc.cores = 40
# )

# ## Human RNA-seq ---------------------------
# 'python3 ./src/tAge.py  "./GEO_hs_rna/processed_data" "./GEO_hs_rna/sample_age" "human"'

# ## Mouse RNA-seq --------------------------
# # Run the tAge program for mouse RNA-seq data.
# # The first argument is the path to the processed data.
# # The second argument is the path to the sample age information.
# # The third argument is the species name.
# 'python3 ./src/tAge.py "./GEO_mm_rna/processed_data" "./GEO_mm_rna/sample_age" "mouse"'
