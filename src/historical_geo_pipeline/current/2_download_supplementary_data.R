library(pbmcapply)
## This step is only required for the RNA-seq data and mouse DNAm data
options(timeout = Inf)
download_and_untar_gse_supp = function(gseid, folder) {
  target_dir = glue("{folder}/supplementary_file_untar/{gseid}")
  if (
    !file.exists(glue("{folder}/supplementary_file/{gseid}_RAW.tar")) || !file.exists(target_dir)
  ) {
    message(glue("Downloading {gseid}"))
    geograbi::geograbi.download.supplementary.file(
      glue("{folder}/supplementary_file"),
      gseid,
      glue("{gseid}_RAW.tar")
    )
  } else {
    message(glue("{gseid} already downloaded"))
  }
  if (! file.exists(target_dir) && file.exists(glue("{folder}/supplementary_file/{gseid}_RAW.tar"))) {
    dir.create(target_dir, recursive = T)
    untar(
      glue("{folder}/supplementary_file/{gseid}_RAW.tar"),
      exdir = target_dir
    )
  }
}
download_and_untar_gse_supp = possibly(download_and_untar_gse_supp, otherwise = NULL)

gse_rna_hs = read_csv("gse_rna_hs.csv")
gse_rna_mm = read_csv("gse_rna_mm.csv")
# gse_me_mm = read_csv("gse_me_mm.csv")

# hs_rna
all_accession = gse_rna_hs$accession

pbmclapply(
  all_accession,
  download_and_untar_gse_supp,
  folder = "GEO_hs_rna",
  mc.cores = 70
)
# GSE97693

# mouse_rna
all_accession = gse_rna_mm$accession

pbmclapply(
  all_accession,
  download_and_untar_gse_supp,
  folder = "GEO_mm_rna",
  mc.cores = 50
)
# GSE15

# mouse_dnam

# files = dir("../02_GEO_DNAm/GEO_cache_mouse", full.names = T)
# files_from = files[str_detect(files, "GSE[0-9]*.tar")]
# pbmclapply(
#   files_from,
#   function(x) {
#     message(x)
#     file_dest = glue("./GEO_mm_dnam/supplementary_file/{basename(x) |> str_replace('.tar', '_RAW.tar')}")
#     file.copy(x, file_dest, overwrite = T)
#   },
#   mc.cores = 50
# )

# all_accession = gse_me_mm$accession

# pbmclapply(
#   all_accession,
#   download_and_untar_gse_supp,
#   folder = "GEO_mm_dnam",
#   mc.cores = 50
# )

# Check the downloaded files
# hs_rna
message("hs_rna")
n_all_accession = length(gse_rna_hs$accession)
all_download = dir("GEO_hs_rna/supplementary_file", full.names = T)
all_attempt = dir("GEO_hs_rna/supplementary_file_untar", full.names = T)
message("N accessions:", n_all_accession)
message("Downloaded: ", length(all_download))
message("Attepts: ", length(all_attempt))

# mouse_rna
message("mouse_rna")
n_all_accession = length(gse_rna_mm$accession)
all_download = dir("GEO_mm_rna/supplementary_file", full.names = T)
all_attempt = dir("GEO_mm_rna/supplementary_file_untar", full.names = T)
message("N accessions:", n_all_accession)
message("Downloaded: ", length(all_download))
message("Attepts: ", length(all_attempt))

# mouse_dnam
message("mouse_dnam")
n_all_accession = length(gse_me_mm$accession)
all_download = dir("GEO_mm_dnam/supplementary_file", full.names = T)
all_attempt = dir("GEO_mm_dnam/supplementary_file_untar", full.names = T)
message("N accessions:", n_all_accession)
message("Downloaded: ", length(all_download))
message("Attepts: ", length(all_attempt))

raw_meta_files = dir("GEO_hs_rna/metadata", full.names = T)
raw_meta_file = raw_meta_files[999]
raw_meta = read_csv(raw_meta_file)
