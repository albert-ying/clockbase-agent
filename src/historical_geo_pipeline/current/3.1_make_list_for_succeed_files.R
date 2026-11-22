tidylog()

# human RNA-seq

file1 = dir("./GEO_hs_rna/archs4/human_data/")
file2 = dir("./GEO_hs_rna/processed_data/", pattern = "csv")

all_meta = read_csv("geo.csv")

success1 = file1 |>
  str_extract("GSE[0-9]+") |>
  unique()
success2 = file2 |>
  str_extract("GSE[0-9]+") |>
  unique()

has_meta = all_meta$accession |> unique()
(success1 %in% has_meta) |> sum()
(!success1 %in% has_meta) |> sum()
(success2 %in% has_meta) |> sum()
(!success2 %in% has_meta) |> sum()

success1[!success1 %in% has_meta][1]

success1

all_success_human = c(success1, success2) |>
  unique()

all_gse = read_csv("gse_rna_hs.csv")
all_gse |>
  # filter(accession %in% all_success_human) |>
  filter(accession %in% success1) |>
  select(accession, samples) |>
  summarise(n = sum(samples, na.rm = T))
write_csv(tibble(id = all_success_human), "all_success_human_RNA.csv")

# mouse RNA-seq

file1 = dir("./GEO_mm_rna/archs4/mouse_data/")
file2 = dir("./GEO_mm_rna/processed_data/", pattern = "csv")

success1 = file1 |>
  str_extract("GSE[0-9]+") |>
  unique()
success2 = file2 |>
  str_extract("GSE[0-9]+") |>
  unique()

all_success_mouse = c(success1, success2) |>
  unique()

all_gse = read_csv("gse_rna_mm.csv")
all_gse |>
  # filter(accession %in% all_success_human) |>
  filter(accession %in% success1) |>
  select(accession, samples) |>
  summarise(n = sum(samples, na.rm = T))

write_csv(tibble(id = all_success_mouse), "all_success_mouse_RNA.csv")
