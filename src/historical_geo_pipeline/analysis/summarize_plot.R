library(parapurrr)
library(ggrepel)
tidylog()

human_dnam = read_csv("result/qced_human_dnam_age.csv")
human_rna = read_csv("result/human_rna_age.csv")
mouse_rna = read_csv("result/mouse_rna_age.csv")

human_dnam_geo = unique(human_dnam$gse_id)
mouse_rna_geo = unique(mouse_rna$gse_id)
human_rna_geo = unique(human_rna$gse_id)

## number of dataset
# Verbally print the amount of datasets, use paste0
print(paste0("Human DNAm: ", length(human_dnam_geo)))
print(paste0("Total DNAm: ", length(mouse_dnam_geo) + length(human_dnam_geo)))
print(paste0("Mouse RNA: ", length(mouse_rna_geo)))
print(paste0("Human RNA: ", length(human_rna_geo)))
print(paste0("Total RNA: ", length(mouse_rna_geo) + length(human_rna_geo)))
print(
  paste0(
    "Total: ",
     length(human_dnam_geo) + length(mouse_rna_geo) + length(human_rna_geo)
  ))

# [1] "Mouse DNAm: 38"

## Make N dataset table ------

geom_donut = function(ds_tb) {
  ## First two clumns are "category" and "value"
  colnames(ds_tb)[1:2] = c("cate", "n")
  ds_tb$fraction <- ds_tb$n / sum(ds_tb$n)
  ds_tb$ymax <- cumsum(ds_tb$fraction)
  ds_tb$ymin <- c(0, head(ds_tb$ymax, n=-1))
  ds_tb$labelPosition <- (ds_tb$ymax + ds_tb$ymin) / 2
  ds_tb$label <- paste0(ds_tb$cate, "\n", format(ds_tb$n, big.mark = ",", scientific = FALSE))

  ggplot(ds_tb, aes(ymax=ymax, ymin=ymin, xmax=2, xmin=0.5, fill=cate, group = species)) +
  # ggplot(ds_tb, aes(ymax=ymax, ymin=ymin, xmax=c(2, 2, 4, 4), xmin=c(0.5, 0.5, 2.5, 2.5), fill=cate)) +
    geom_rect(fill = "black", xmin = 0.4, xmax = 2.0) +
    geom_rect() +
    # geom_label(x=4, aes(y=labelPosition, label=label, fill=cate), size=6, alpha = 0.8) + # x here controls label position (inner / outer)
    geom_text(x=3, aes(y=labelPosition, label=label, color=cate), size=5) + # x here controls label position (inner / outer)
    coord_polar(theta="y") +
    xlim(c(-1, 4)) +
    theme_void() +
    theme(legend.position = "none")
}

ds_tb = tibble(
  species = c("human", "mouse", "human"),
  type = c("DNAm", "RNA", "RNA"),
  n = c(length(human_dnam_geo), length(mouse_rna_geo), length(human_rna_geo))
) |>
  mutate(cate = glue("{type} {species}")) |>
  select(cate, n, everything())

p1 = ds_tb |>
  geom_donut()
  # scale_fill_paletteer_d("nord::frost") +
  # scale_color_paletteer_d("nord::frost")

ggsave('plot/dataset_summary.pdf', p1, width = 8, height = 8)
 
## number of samples
mouse_dnam = 1749
human_dnam = nrow(human_dnam)
mouse_rna = nrow(mouse_rna)
human_rna = nrow(human_rna)
 
## verbally print the amount of samples
print(paste0("Mouse DNAm: ", mouse_dnam))
print(paste0("Human DNAm: ", human_dnam))
print(paste0("Total DNAm: ", mouse_dnam + human_dnam))
print(paste0("Mouse RNA: ", mouse_rna))
print(paste0("Human RNA: ", human_rna))
print(paste0("Total RNA: ", mouse_rna + human_rna))
print(paste0("Total: ", mouse_dnam + human_dnam + mouse_rna + human_rna))
 
# [1] "Mouse DNAm: 1749"
# [1] "Human DNAm: 192635"
# [1] "Total DNAm: 194384"
# [1] "Mouse RNA: 47961"
# [1] "Human RNA: 30346"
# [1] "Total RNA: 78307"
# [1] "Total: 272691"

sample_tb = tibble(
  species = c("mouse", "human", "mouse", "human"),
  type = c("DNAm", "DNAm", "RNA", "RNA"),
  n = c(mouse_dnam, human_dnam, mouse_rna, human_rna)
) |>
  mutate(cate = glue("{type} {species}")) |>
  select(cate, n, everything())

p2 = sample_tb |>
  geom_donut()

ggsave('plot/sample_summary.pdf', p2, width = 8, height = 8)
