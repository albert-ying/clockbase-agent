library(ClockBasis)
library(geograbi)
tidylog()

## Human DNAm --------------------
## Make clock list ---------------

files = c(
  "Elasticnet_alltissues_5K_coefs_fixed.csv",
  "Elasticnet_blood_5K_coefs_fixed.csv"
)

path = "./mortality_clocks/"

clock_ls = str_c(path, files) |>
  map(read_csv) |>
  map(~filter(.x, coefs != 0)) |>
  set_names(c("MortClock_multi", "MortClock_blood"))

run_human_clock = function(file, out_folder) {
  gse_id = str_extract(file, "GSE[0-9]+")

  outfile = glue("{out_folder}/{gse_id}.csv")

  if (file.exists(outfile)) {
    return(NULL)
  }

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
    ClockBasis::do_dnam_clock_human(clock = clock_ls, ncores = 7)

  write_csv(dat, outfile)
}

possible_run_human_clock = possibly(run_human_clock, otherwise = NULL, quiet = FALSE)

if (sys.nframe() == 0) {
  files = dir("./GEO_hs_dnam/series_matrix/", full.names = T)
  out_folder = "./GEO_hs_dnam/sample_age_mort"
  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }

  pbmcapply::pbmclapply(
    files,
    possible_run_human_clock,
    out_folder = out_folder,
    mc.cores = 15
  )
}
