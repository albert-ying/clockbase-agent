library(ClockBasis)
library(geograbi)

## Human DNAm --------------------


## Make clock list ---------------

causal_files = c(
  "Aging-GIP1-adj_cau_eln_clock_padjbonf_rank_weight_force1_03",
  "Aging-GIP1-adj_dam_eln_clock_padjbonf_rank_weight_force1_03",
  "Aging-GIP1-adj_pro_eln_clock_padjbonf_rank_weight_force1_03"
)

causal_path = "/Lasagna/kying/fasrc/1_causal_aging_omic/new_clock_models5/"

causal_clock_ls = str_c(causal_path, causal_files, ".csv") |>
  map(read_csv) |>
  set_names(c("CausAge", "DamAge", "AdaptAge"))

causal_clock = \(data, ...) {
  ClockBasis::do_dnam_clock_human(data, clock = causal_clock_ls, ...)
}

run_human_clock = function(file, out_folder) {
  gse_id = str_extract(file, "GSE[0-9]+")

  mat = geograbi::geograbi.read.gse.matrix(file)

  dat = mat$data |>
    as.data.frame() |>
    rownames_to_column("probe_id") |>
    as_tibble() |>
    ClockBasis::do_dnam_clock_human(clock = causal_clock_ls)

  write_csv(dat, glue("{out_folder}/{gse_id}.csv"))
}

possible_run_human_clock = possibly(run_human_clock, otherwise = NULL)

if (sys.nframe() == 0) {
  files = dir("./GEO_hs_dnam/series_matrix/", full.names = T)
  out_folder = "./GEO_hs_dnam/sample_age_causal"
  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }

  pbmcapply::pbmclapply(
    files,
    possible_run_human_clock,
    out_folder = out_folder,
    mc.cores = 10
  )
}
