library(ggrepel)
library(pbmcapply)
tidylog()

#- Predefined functions ----
test_plot = function(dem_sum, group_var, value_var, fill_var) {
  require(rstatix)
  require(ggpubr)
  require(gridExtra)
  require(grid)
  require(patchwork)
  add_pval = function(data) {
    tryCatch(
      {
        stat_pvalue_manual(
          data,
          label = "p.adj.signif",
          inherit.aes = F,
          hide.ns = T,
          label.size = 8,
          bracket.size = 0.6,
          tip.length = 0.02,
          step.increase = 0.05,
          bracket.nudge.y = 1,
          vjust = 0.65 
        )
      },
      error = function(err) {
     }
    )
  }
  df = dem_sum[, c(group_var, value_var, fill_var)]
  colnames(df) = c("x", "y", "c")
  df = df |>
    mutate(x = as.factor(as.character(x)))
  res = df |>
    anova_test(y ~ x) |>
    as_tibble()
  # Perform T test
  res_t = df |>
    pairwise_t_test(y ~ x) |>
    adjust_pvalue(method = "holm") |>
    add_significance() |>
    add_y_position(scales = "free", step.increase = 0) |>
    filter(p<0.05)
  # Ploting
  df = df |>
    cbind(res) |>
    mutate(anova = glue("ANOVA P = {p}")) |>
    mutate(anova = ifelse(p < 0.05, glue("<b style='color:orange'>{anova}</b>"), anova)) 
  p = {df |>
    ggplot(aes(x = x, y = y, fill = c)) +
    geom_violin(scale = "width", width = 0.5, alpha = 0.3, position=position_dodge(width = 0.6), trim = FALSE) +
    labs(x = group_var, y = value_var, fill = fill_var, caption = unique(df$anova))
    } |>
    base_mode() +
    geom_boxplot(width = 0.1, alpha = 0.5, position=position_dodge(width = 0.6)) +
    geom_sina(alpha = 0.8, pch = 21, position=position_dodge(width = 0.6)) +
    add_pval(res_t)

  if (fill_var == group_var) {
    p = p + theme(legend.position = "none")
  }

  return(p)
}

clocks = c(
  "Horvathv1", "Hannum", "Lin", "PhenoAge", "YingCausAge", "YingDamAge",
  "YingAdaptAge", "Horvathv2", "PEDBE", "Zhang_10", "DunedinPoAm38", "DunedinPACE",
  "GrimAgeV1", "GrimAgeV2", "DNAmTL", "Knight", "LeeRobust"
)

biomarkers = c(
  "DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS",
  "DNAmlogA1C", "DNAmlogCRP", "Alcohol", "BMI", "Smoking"
)

blood_clocks = c(
  "Horvathv1", "Hannum", "Lin", "PhenoAge", "YingCausAge", "YingDamAge",
  "YingAdaptAge", "Horvathv2", "GrimAgeV1", "GrimAgeV2"
)

human_dnam = read_csv("result/qced_human_dnam_age.csv")

human_dnam_meta = read_csv("result/gpt4_human_dnam_metadata.csv")


blood_samples = human_dnam_meta |>
  filter(tissue == "Whole Blood") |>
  pull(gsm_id)

human_dnam_cleaned = human_dnam |>
  pivot_longer(cols = -c(gse_id, gsm_id), names_to = "clock", values_to = "DNAm_Age") |>
  group_by(clock) |>
  drop_na(DNAm_Age) |>
  filter(!abs(DNAm_Age - median(DNAm_Age, na.rm = TRUE)) > 1.5 * IQR(DNAm_Age, na.rm = TRUE)) |>
  ungroup() |>
  pivot_wider(names_from = clock, values_from = DNAm_Age)

human_dnam_long = human_dnam_cleaned |>
  pivot_longer(cols = -c(gse_id, gsm_id), names_to = "clock", values_to = "DNAm_Age")

## - Distribution ----
# Function to create density plot
create_density_plot <- function(data, sample_ids, title, filename) {
  p = data |>
    filter(gsm_id %in% sample_ids) |>
    mutate(clock = as.factor(clock)) |>
    ggplot(aes(x = DNAm_Age)) +
    geom_density(aes(fill = clock), alpha = 0.5) +
    labs(x = "Age", y = "Density", title = title) +
    facet_wrap(~ clock, scales = "free") +
    theme(legend.position = "none",
          axis.text.y = element_blank()) +
    paletteer::scale_fill_paletteer_d("pals::glasbey") +
    scale_x_continuous(labels = function(x) {
      ifelse(x > 10000, format(x, scientific = TRUE, drop0trailing = T), x)
    })
  
  ggsave(filename, p, w = 28, h = 10)
}

# Create density plots for all samples and blood samples
create_density_plot(human_dnam_long, unique(human_dnam_long$gsm_id), "All Samples", "plot/human_dnam_age_dist_all.pdf")
create_density_plot(human_dnam_long, blood_samples, "Blood Samples", "plot/human_dnam_age_dist_blood.pdf")

# Function to create ridges plot
create_ridges_plot <- function(data, sample_ids, title, filename) {
  library(ggridges)
  
  p = data |>
    filter(gsm_id %in% sample_ids) |>
    filter(clock %in% blood_clocks) |>
    group_by(clock) |>
    filter(DNAm_Age >= 0 & DNAm_Age <= 100) |>
    mutate(mean_DNAm_Age = median(DNAm_Age, na.rm = T)) |>
    mutate(clock = fct_reorder(clock, mean_DNAm_Age)) |>
    ggplot(aes(x = DNAm_Age, y = clock, fill = after_stat(x))) +
    geom_density_ridges_gradient(alpha = 0.5) +
    scale_fill_viridis_c(name = "Predicted Age", option = "C") +
    labs(y = NULL, title = title)
  
  ggsave(filename, p, w = 10, h = 5)
}

# Create ridges plots for all samples and blood samples
create_ridges_plot(human_dnam_long, unique(human_dnam_long$gsm_id), "All Samples", "plot/human_dnam_age_dist_all_ridges.pdf")
create_ridges_plot(human_dnam_long, blood_samples, "Blood Samples", "plot/human_dnam_age_dist_blood_ridges.pdf")

# - Correlation matrix ------
library(rstatix)
library(corrr)
library(corrplot)
library(purrr)
library(dplyr)

print(n_distinct(human_dnam_cleaned$gse_id))
print(length(human_dnam_cleaned$gsm_id))
# 100157 samples
# 1923 datasets

# Function to compute correlation matrix
compute_cor_matrix <- function(data) {
  data %>%
    correlate(diagonal = 1) %>%
    as.data.frame() %>%
    column_to_rownames("term") %>%
    as.matrix()
}

# Function to compute p-value matrix
compute_p_matrix <- function(data) {
  data %>%
    select(-gse_id, -gsm_id) %>%
    colpair_map(function(a, b) {
      cor.test(a, b)$p.value
    }, .diagonal = 1) %>%
    as.data.frame() %>%
    column_to_rownames("term") %>%
    as.matrix()
}

# Function to add prefix to biomarkers
add_prefix_to_biomarkers <- function(matrix, biomarkers, prefix = "*") {
  for (i in seq_len(nrow(matrix))) {
    if (rownames(matrix)[i] %in% biomarkers) {
      rownames(matrix)[i] <- paste0(prefix, rownames(matrix)[i])
      colnames(matrix)[i] <- paste0(prefix, colnames(matrix)[i])
    }
  }
  matrix
}

# Function to plot correlation matrix
plot_cor_matrix <- function(cor_mat, p_mat, pthres, filename) {
  pdf(filename, width = 11, height = 7)
  corrplot(
    cor_mat,
    method = "square",
    outline = TRUE,
    order = "hclust",
    addgrid.col = "darkgray",
    type = "full",
    tl.col = "black",
    tl.pos = "tp",
    diag = TRUE,
    p.mat = p_mat,
    insig = "label_sig",
    sig.level = pthres,
    pch = "·",
    pch.cex = 1,
    pch.col = "black",
    na.label = "square",
    na.label.col = "white"
  )
  dev.off()
}
# Function to process and plot correlation matrices
process_and_plot_cor_matrices <- function(data, suffix = "") {
  cor_mat <- compute_cor_matrix(data)
  p_mat <- compute_p_matrix(data)
  
  write.csv(cor_mat, paste0("result/human_dnam_cor_mat", suffix, ".csv"))
  write.csv(p_mat, paste0("result/human_dnam_p_mat", suffix, ".csv"))
  
  pthres <- 0.05 / ((dim(p_mat)[1] * dim(p_mat)[2] - dim(p_mat)[1]) / 2)
  
  cor_mat <- add_prefix_to_biomarkers(cor_mat, biomarkers)
  p_mat <- add_prefix_to_biomarkers(p_mat, biomarkers)

  plot_cor_matrix(cor_mat, p_mat, pthres, paste0("plot/human_corr_mat", suffix, ".pdf"))
}

# Process all samples
process_and_plot_cor_matrices(human_dnam_cleaned)

# Process blood samples only
blood_samples_df <- human_dnam_cleaned %>%
  filter(gsm_id %in% blood_samples)

print(paste("Number of unique datasets:", n_distinct(blood_samples$gse_id)))
print(paste("Number of samples:", length(blood_samples$gsm_id)))

process_and_plot_cor_matrices(blood_samples_df, "_blood")

## - PCA for clock cluster ----
mat = human_dnam_cleaned |>
  # select(-gse_id, -gsm_id) |>
  select(all_of(clocks)) |>
  drop_na() |>
  as.matrix() |>
  t()

pca_res = prcomp(mat, scale. = TRUE, center = TRUE)

var_explained = pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

pp = pca_res |>
  broom::tidy() |>
  filter(PC %in% 1:2) |>
  pivot_wider(names_from = PC, values_from = value) |>
  ggplot(aes(`1`, `2`)) |>
  geom_point(aes(fill = row), show.legend = F) |>
  geom_text_repel(aes(label = row), max.overlaps = Inf, size = 5) |>
  scale_fill_paletteer_d("pals::stepped") |>
  base_mode() +
  labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]))
  
ggsave2("plot/clock_pca.pdf", pp, w = 8, h = 6)

# PCA for blood samples only ----
blood_mat = blood_samples_df |>
  select(all_of(clocks)) |>
  drop_na() |>
  as.matrix() |>
  t()

blood_pca_res = prcomp(blood_mat, scale. = TRUE, center = TRUE)

blood_var_explained = blood_pca_res$sdev^2 / sum(blood_pca_res$sdev^2) * 100

blood_pp = blood_pca_res |>
  broom::tidy() |>
  filter(PC %in% 1:2) |>
  pivot_wider(names_from = PC, values_from = value) |>
  ggplot(aes(`1`, `2`)) |>
  geom_point(aes(fill = row), show.legend = F) |>
  geom_text_repel(aes(label = row), max.overlaps = Inf, size = 5) |>
  scale_fill_paletteer_d("pals::stepped") |>
  base_mode() +
  labs(x = sprintf("PC1 (%.1f%%)", blood_var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", blood_var_explained[2]))
  
ggsave2("plot/human_blood_clock_pca.pdf", blood_pp, w = 8, h = 6)



##- Sample PCA and UMAP ---------
data = recipes::recipe(formula(~.), data = select(human_dnam_cleaned, gse_id, gsm_id, all_of(clocks))) %>%
  recipes::step_impute_mean(recipes::all_numeric_predictors()) %>%
  recipes::step_normalize(recipes::all_numeric_predictors()) %>%
  embed::step_umap(recipes::all_numeric_predictors(), neighbors = 200, num_comp = 3) %>%
  recipes::prep() %>%
  recipes::juice()

write_csv(data, "result/human_dnam_umap_age.csv")

# data = read_csv("result/human_dnam_umap_age.csv")
plot_data = data |>
  inner_join(human_dnam_cleaned)

annotation = tribble(
  ~gsm_id, ~label,
  "GSM1589932", "hESC/iPSC",
  "GSM1854190", "HepaRG differentiation",
  "GSM1439952", "Infant blood",
  "GSM1842947", "Placenta",
  "GSM4080508", "Pediatric buccal cell",
  "GSM2467060", "Fetal brain",
  "GSE110184", "Breast cancer",
  "GSM643345", "DNMT over-expressed cell",
  "GSM1873846", "Esophageal adenocarcinoma",
  "GSM1818245", "Germ cell tumor",
  "GSM1342235", "Pediatric glioblastomas",
  "GSM4151384", "Sperm",
  "GSM1505332", "Nuclear-transfer ESC",
  "GSM2256793", "Cord blood"
) |>
  inner_join(plot_data) |>
  select(gsm_id, label, UMAP1, UMAP2) |>
  distinct()


# get the rounded 75% range
plotlist = list()
clocks_to_plot = c(
  "Horvathv1", "PhenoAge", "YingDamAge",
  "YingAdaptAge", "GrimAgeV2", "DunedinPACE"
)

for (i in 1:length(clocks_to_plot)) {
  clock = clocks_to_plot[i]

  quants = plot_data[[clock]] |>
    quantile(probs = c(0.05, 0.95), na.rm = T)

  if (quants[[2]] > 10) {
    q1 = round(quants[[1]] / 10) * 10
    q2 = round(quants[[2]] / 10) * 10
    if (q1 == q2) {
      q1 = q1 - 5
      q2 = q2 + 5
    }
  } else {
    q1 = round(quants[[1]] * 2) / 2
    q2 = round(quants[[2]] * 2) / 2
  }

  plot_data |>
    filter(UMAP2 > 4) |>
    select(gsm_id) |>
    pull()

  plot = plot_data |>
    arrange(UMAP1, UMAP2) |>
    # remove outliers in the plot
    filter(!abs(UMAP1 - median(UMAP1)) > 3 * sd(UMAP1) & !abs(UMAP2 - median(UMAP2)) > 3 * sd(UMAP2)) |>
    mutate(DNAmAge = ifelse(!!sym(clock) > q2, q2, !!sym(clock))) |>
    mutate(DNAmAge = ifelse(DNAmAge < q1, q1, DNAmAge)) |>
    ggplot(aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = DNAmAge), shape = 19, alpha = 1, size = 0.1) +
    labs(color = clock) +
    theme_void() +
    theme(legend.position = "top") +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    better_color_legend

  plotlist[[i]] = plot

  plot2 = plot +
    geom_mark_ellipse(
      data = annotation,
      aes(label = label, group = label),
      label.fill = "#ffffffCC",
      alpha = 0.6,
      expand = unit(.1, "mm"),
      con.cap = unit(.0, "mm")
    )
  ggsave(glue("plot/human_umap_{clock}.png"), plot2, w = 15, h = 10, bg = "white", dpi = 600)
}

tot_plot = plotlist |> 
  wrap_plots(nrow = 2) * theme(axis.text.x = element_blank(), axis.text.y = element_blank())
ggsave("plot/human_dnam_umap_panels.png", tot_plot, w = 18, h = 8, bg = "white", dpi = 600)

## Write the plotly obj
library(plotly)
plotly_umap = function(pp) {
  plot_ly(
    x = pp$UMAP1,
    y = pp$UMAP2,
    z = pp$UMAP3,
    type = "scatter3d",
    hoverinfo = "text",
    text = pp$label,
    marker = list(
      color = pp$DNAmAge,
      colorscale = "Portland",
      alpha = 0.6,
      size = 1,
      showscale = TRUE
    )
  ) |>
    plotly::layout(
      title = clock,
      scene = list(
        title = clock,
        xaxis = list(title = "UMAP1"),
        yaxis = list(title = "UMAP2"),
        zaxis = list(title = "UMAP3")
      ),
      paper_bgcolor = "#000",
      plot_bgcolor = "#000",
      font = list(color = "#fff")
    )
}

for (clock in clocks_to_plot) {
  quants = plot_data[[clock]] |>
    quantile(probs = c(0.05, 0.95), na.rm = T)

  if (quants[[2]] > 10) {
    q1 = round(quants[[1]] / 10) * 10
    q2 = round(quants[[2]] / 10) * 10
    if (q1 == q2) {
      q1 = q1 - 5
      q2 = q2 + 5
    }
  } else {
    q1 = round(quants[[1]] * 2) / 2
    q2 = round(quants[[2]] * 2) / 2
  }

  pp = plot_data |>
    arrange(UMAP1, UMAP2, UMAP3) |>
    filter(
      abs(UMAP1 - median(UMAP1)) < 3 * sd(UMAP1) & abs(UMAP2 - median(UMAP2)) < 3 * sd(UMAP2) &
      abs(UMAP3 - median(UMAP3)) < 3 * sd(UMAP3)
    ) |>
    left_join(human_dnam_meta, by = c("gsm_id", "gse_id")) |>
    mutate(annotation_lite = glue(
      "</br>{gse_id} {gsm_id} </br> Age: {age} </br> Sex: {sex} </br> Tissue: {tissue} </br> Disease: {disease}"
    )) |>
    mutate(DNAmAge = ifelse(!!sym(clock) > q2, q2, !!sym(clock))) |>
    mutate(DNAmAge = ifelse(DNAmAge < q1, q1, DNAmAge)) |>
    mutate(label = str_c("</br>Predicted Age: ", signif(!!sym(clock), 3), annotation_lite))

  plotly_plot = plotly_umap(pp)

  # saveRDS(plotly_plot, glue("plot/human_dnam_umap_{clock}.rds"))
  saveRDS(plotly_plot, glue("plotly/human_dnam_umap_{clock}.rds"))
  htmlwidgets::saveWidget(plotly_plot, glue("plotly/human_dnam_umap_{clock}.html"), selfcontained = T)
}

# Blood only version
blood_data = recipes::recipe(formula(~.), data = select(blood_samples_df, gse_id, gsm_id, all_of(clocks))) %>%
  recipes::step_impute_mean(recipes::all_numeric_predictors()) %>%
  recipes::step_normalize(recipes::all_numeric_predictors()) %>%
  embed::step_umap(recipes::all_numeric_predictors(), neighbors = 200, num_comp = 3) %>%
  recipes::prep() %>%
  recipes::juice()

write_csv(blood_data, "result/human_blood_dnam_umap_age.csv")

blood_plot_data = blood_data |>
  inner_join(blood_samples_df)

blood_plotlist = list()

for (i in 1:length(clocks_to_plot)) {
  clock = clocks_to_plot[i]

  quants = blood_plot_data[[clock]] |>
    quantile(probs = c(0.05, 0.95), na.rm = T)

  if (quants[[2]] > 10) {
    q1 = round(quants[[1]] / 10) * 10
    q2 = round(quants[[2]] / 10) * 10
    if (q1 == q2) {
      q1 = q1 - 5
      q2 = q2 + 5
    }
  } else {
    q1 = round(quants[[1]] * 2) / 2
    q2 = round(quants[[2]] * 2) / 2
  }

  blood_plot = blood_plot_data |>
    arrange(UMAP1, UMAP2) |>
    filter(!abs(UMAP1 - median(UMAP1)) > 3 * sd(UMAP1) & !abs(UMAP2 - median(UMAP2)) > 3 * sd(UMAP2)) |>
    mutate(DNAmAge = ifelse(!!sym(clock) > q2, q2, !!sym(clock))) |>
    mutate(DNAmAge = ifelse(DNAmAge < q1, q1, DNAmAge)) |>
    ggplot(aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = DNAmAge), shape = 19, alpha = 1, size = 0.1) +
    labs(color = clock) +
    theme_void() +
    theme(legend.position = "top") +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    better_color_legend

  blood_plotlist[[i]] = blood_plot

  ggsave(glue("plot/human_blood_umap_{clock}.png"), blood_plot, w = 15, h = 10, bg = "white", dpi = 600)
}

blood_tot_plot = blood_plotlist |> 
  wrap_plots(nrow = 2) * theme(axis.text.x = element_blank(), axis.text.y = element_blank())
ggsave("plot/human_blood_dnam_umap_panels.png", blood_tot_plot, w = 18, h = 8, bg = "white", dpi = 600)

for (clock in clocks_to_plot) {
  quants = blood_plot_data[[clock]] |>
    quantile(probs = c(0.05, 0.95), na.rm = T)

  if (quants[[2]] > 10) {
    q1 = round(quants[[1]] / 10) * 10
    q2 = round(quants[[2]] / 10) * 10
    if (q1 == q2) {
      q1 = q1 - 5
      q2 = q2 + 5
    }
  } else {
    q1 = round(quants[[1]] * 2) / 2
    q2 = round(quants[[2]] * 2) / 2
  }

  blood_pp = blood_plot_data |>
    arrange(UMAP1, UMAP2, UMAP3) |>
    filter(
      abs(UMAP1 - median(UMAP1)) < 3 * sd(UMAP1) & abs(UMAP2 - median(UMAP2)) < 3 * sd(UMAP2) &
      abs(UMAP3 - median(UMAP3)) < 3 * sd(UMAP3)
    ) |>
    left_join(human_dnam_meta, by = c("gsm_id", "gse_id")) |>
    mutate(annotation_lite = glue(
      "</br>{gse_id} {gsm_id} </br> Age: {age} </br> Sex: {sex} </br> Disease: {disease}"
    )) |>
    mutate(DNAmAge = ifelse(!!sym(clock) > q2, q2, !!sym(clock))) |>
    mutate(DNAmAge = ifelse(DNAmAge < q1, q1, DNAmAge)) |>
    mutate(label = str_c("</br>Predicted Age: ", signif(!!sym(clock), 3), annotation_lite))

  blood_plotly_plot = plotly_umap(blood_pp)

  saveRDS(blood_plotly_plot, glue("plotly/human_blood_dnam_umap_{clock}.rds"))
  htmlwidgets::saveWidget(blood_plotly_plot, glue("plotly/human_blood_dnam_umap_{clock}.html"), selfcontained = T)
}

## =================================================

## Analyze the age Acceleration! (with regression) --------------------

# no_accelerat_clocks = c("DunedinPoAm", "DunedinPACE")

# 230506 samples
# human_dnam$gse_id |> n_distinct()
# 2653 datasets

age_tb = human_dnam_meta |>
  select(1:2, harmonized.age) |>
  filter(harmonized.age >= 0 & harmonized.age <= 120)
# 80,346 samples
# age_tb$GEO |> n_distinct()
# 618 datasets

accer_tb = human_dnam_long |>
  separate(sample_id, c("GEO", "GSM"), sep = "_") |>
  inner_join(age_tb, by = c("GSM", "GEO")) |>
  filter(DNAm_Age_missing <= 0.2)
  # drop_na(harmonized.age) |>
  # mutate(value = ifelse(name %in% no_accelerat_clocks, value, value - harmonized.age))

accer_tb = accer_tb |>
  group_by(clock) |>
  group_map(~ {
    res = resid(lm(DNAm_Age ~ harmonized.age, data = .x))
    .x$accel = res
    .x
  }, .keep = T) |>
  bind_rows()

human_dnam_accer_cl = accer_tb |>
  select(sample_id = GSM, clock, accel) |>
  distinct() |>
  pivot_wider(id_cols = sample_id, names_from = clock, values_from = accel, values_fn = mean)

## Correlation matrix ------
library(rstatix)
library(corrr)
library(corrplot)

cor_mat = human_dnam_accer_cl |>
  correlate(diagonal = 1) |>
  as.data.frame() |>
  column_to_rownames("term") |>
  as.matrix()

p_mat = human_dnam_accer_cl |>
  select(-sample_id) |>
  colpair_map(function(a, b) {
    cor.test(a, b)$p.value
  }, .diagonal = 1) |>
  as.data.frame() |>
  column_to_rownames("term") |>
  as.matrix()

write_csv(as.data.frame(cor_mat), "data/human_dnam/cor_mat_resid.csv")
write_csv(as.data.frame(p_mat), "data/human_dnam/p_mat_resid.csv")

pthres = 0.05 / ((dim(p_mat)[1] * dim(p_mat)[2] - dim(p_mat)[1]) / 2)

# colnames(cor_mat) = ifelse(
#   colnames(cor_mat) %in% no_accelerat_clocks,
#   colnames(cor_mat),
#   paste(colnames(cor_mat), "-accel")
# )
# rownames(cor_mat) = ifelse(
#   rownames(cor_mat) %in% no_accelerat_clocks,
#   rownames(cor_mat),
#   paste(rownames(cor_mat), "-accel")
# )
# colnames(p_mat) = colnames(cor_mat)
# rownames(p_mat) = rownames(cor_mat)

corrplot(cor_mat)

pdf("plot/human_corr_mat_accel.pdf", width = 12, height = 8)
corrplot(
  cor_mat,
  method = "square",
  outline = T,
  order = "hclust",
  addgrid.col = "darkgray",
  type = "full",
  tl.col = "black",
  tl.pos = "tp",
  diag = T,
  p.mat = p_mat,
  insig = "label_sig",
  sig.level = pthres,
  pch = "·",
  pch.cex = 1,
  pch.col = "black",
  na.label = "square",
  na.label.col = "white"
)
# corrplot(
#   cor_mat,
#   method = "number",
#   outline = T,
#   order = "hclust",
#   addgrid.col = "darkgray",
#   type = "lower",
#   diag = F,
#   col = "black",
#   number.cex = 0.5,
#   tl.pos = "n",
#   cl.pos = "n",
#   pch.cex = 1.5,
#   pch.col = "white",
#   na.label = "square",
#   na.label.col = "white",
#   add = T
# )
dev.off()

## Analyze the age Acceleration (substraction)! --------------------

no_accelerat_clocks = c("DunedinPoAm", "DunedinPACE")
# age_clock

human_dnam_meta_all = read_csv("data/human_dnam/processed_metadata_all.csv")
# 192,635 samples
# human_dnam_all$GEO |> n_distinct()
# 2416 datasets
human_dnam_long = read_csv("data/human_dnam/sample_age_long.csv")

age_tb = human_dnam_meta_all |>
  select(1:2, harmonized.age) |>
  filter(harmonized.age >= 0 & harmonized.age <= 120)
# 80,346 samples
# age_tb$GEO |> n_distinct()
# 618 datasets

accer_tb = human_dnam_long |>
  filter(clock %in% c(age_clock, no_accelerat_clocks)) |>
  separate(sample_id, c("GEO", "GSM"), sep = "_") |>
  inner_join(age_tb, by = c("GSM", "GEO")) |>
  filter(DNAm_Age_missing <= 0.2) |>
  mutate(accel = ifelse(clock %in% age_clock, DNAm_Age - harmonized.age, DNAm_Age))

human_dnam_accer_cl = accer_tb |>
  select(sample_id = GSM, clock, accel) |>
  distinct() |>
  pivot_wider(id_cols = sample_id, names_from = clock, values_from = accel, values_fn = mean)

## Correlation matrix ------
library(rstatix)
library(corrr)
library(corrplot)

cor_mat = human_dnam_accer_cl |>
  correlate(diagonal = 1) |>
  as.data.frame() |>
  column_to_rownames("term") |>
  as.matrix()

p_mat = human_dnam_accer_cl |>
  select(-sample_id) |>
  colpair_map(function(a, b) {
    cor.test(a, b)$p.value
  }, .diagonal = 1) |>
  as.data.frame() |>
  column_to_rownames("term") |>
  as.matrix()
write_csv(as.data.frame(cor_mat), "data/human_dnam/cor_mat_accel2.csv")
write_csv(as.data.frame(p_mat), "data/human_dnam/p_mat_accel2.csv")

pthres = 0.05 / ((dim(p_mat)[1] * dim(p_mat)[2] - dim(p_mat)[1]) / 2)

# colnames(cor_mat) = ifelse(
#   colnames(cor_mat) %in% no_accelerat_clocks,
#   colnames(cor_mat),
#   paste(colnames(cor_mat), "-accel")
# )
# rownames(cor_mat) = ifelse(
#   rownames(cor_mat) %in% no_accelerat_clocks,
#   rownames(cor_mat),
#   paste(rownames(cor_mat), "-accel")
# )
# colnames(p_mat) = colnames(cor_mat)
# rownames(p_mat) = rownames(cor_mat)

corrplot(cor_mat)

pdf("plot/human_corr_mat_accel2.pdf", width = 8, height = 6)
corrplot(
  cor_mat,
  method = "square",
  outline = T,
  order = "hclust",
  addgrid.col = "darkgray",
  type = "full",
  tl.col = "black",
  tl.pos = "tp",
  diag = T,
  p.mat = p_mat,
  insig = "label_sig",
  sig.level = pthres,
  pch = "·",
  pch.cex = 1,
  pch.col = "black",
  na.label = "square",
  na.label.col = "white"
)
# corrplot(
#   cor_mat,
#   method = "number",
#   outline = T,
#   order = "hclust",
#   addgrid.col = "darkgray",
#   type = "lower",
#   diag = F,
#   col = "black",
#   number.cex = 0.5,
#   tl.pos = "n",
#   cl.pos = "n",
#   pch.cex = 1.5,
#   pch.col = "white",
#   na.label = "square",
#   na.label.col = "white",
#   add = T
# )
dev.off()

dens_df = accer_tb |>
  mutate(accel = ifelse(DNAm_Age_missing > 0.2, NA, accel)) |>
  mutate(clock = ifelse(clock %in% no_accelerat_clocks, clock, paste(clock, "accel"))) |>
  mutate(ref = ifelse(clock %in% no_accelerat_clocks, 1, 0))

pp = dens_df |>
  mutate(dage = accel - ref) |>
  group_by(clock) |>
  mutate(mean_dage = mean(dage, na.rm = T)) |>
  ungroup() |>
  mutate(group = ifelse(dage > 0, "Accelerated", "Decelerated")) |>
  mutate(group = ifelse(sign(mean_dage) == sign(dage), glue("***{group}***"), group)) |>
  mutate(group_order = ifelse(str_detect(group, "Accelerated"), 1, 2)) |>
  mutate(adage = abs(dage)) |>
  group_by(clock) |>
  group_split() |>
  map(~ {
    .x |>
    mutate(group = as.factor(group)) |>
    mutate(group = fct_reorder(group, group_order)) |>
    test_plot("group", "adage", "group") +
      geom_boxplot(width = 0.1) +
      labs(subtitle = .x$name[1], x = "", y = "Abs Delta Age")
  }) |> wrap_plots(nrow = 2)

ggsave("plot/human_dnam_accelerat_violin.pdf", pp, w = 32, h = 12)

df = dens_df |>
  arrange(clock) |>
  group_by(clock) |>
  mutate(acc = sum(accel > ref)) |>
  mutate(dec = sum(accel < ref)) |>
  select(clock, acc, dec, ref) |>
  distinct()

stat_df = df |>
  pivot_longer(acc:dec, names_to = "type", values_to = "accel") |>
  group_map(~ rstatix::chisq_test(x = .$accel)) |>
  bind_rows() |>
  bind_cols(df) |>
  mutate(ratio = acc / dec) |>
  mutate(text = glue("Accel:Decel = {signif(ratio, 2)}\n p {ifelse(p < 1e-12, '< 1e-12', paste('=', signif(p, 3)))}"))

pp = dens_df |>
  inner_join(stat_df, by = c("clock", "ref")) |>
  mutate(clock = paste(clock, "\n\n", text)) |>
  mutate(clock = as.factor(clock)) |>
  ggplot(aes(x = accel)) +
  geom_density(aes(fill = clock), alpha = 0.5, show.legend = F) +
  geom_vline(aes(xintercept = ref), color = "tomato", linetype = "dashed") +
  labs(x = "Delta Age", y = "Density") +
  facet_wrap(~clock, scales = "free", nrow = 2)

ggsave("plot/human_dnam_age_accel_dist.pdf", pp, w = 45, h = 10)

## PCA and UMAP ---------
set.seed(123)
data = recipes::recipe(formula(~.), data = human_dnam_accer_cl) %>%
  recipes::step_impute_knn(recipes::all_numeric_predictors(), neighbors = 5) %>%
  recipes::step_normalize(recipes::all_numeric_predictors()) %>%
  recipes::prep() %>%
  recipes::juice()

tsne_res = Rtsne::Rtsne(
  as.matrix(data[,2:ncol(data)]),
  dims = 3,
  check_duplicates = F,
  perplexity = 30,
  max_iter = 1000
)

data = tsne_res$Y |>
  as_tibble() |>
  set_names(c("tsne1", "tsne2", "tsne3")) |>
  bind_cols(data[,1])

write_csv(data, "data/human_dnam/tsne_age_accel.csv")

plot_data = data |>
  inner_join(human_dnam_accer_cl) |>
  drop_na()

# get the rounded 75% range
i = 0
plotlist = list()
for (clock in unique(colnames(human_dnam_accer_cl)[2:ncol(human_dnam_accer_cl)])) {
  i = i + 1
  quants = plot_data[[clock]] |>
    quantile(probs = c(0.05, 0.95), na.rm = T)

  if (quants[[2]] > 10) {
    q1 = round(quants[[1]] / 10) * 10
    q2 = round(quants[[2]] / 10) * 10
  } else {
    q1 = round(quants[[1]] * 2) / 2
    q2 = round(quants[[2]] * 2) / 2
  }

  plot = plot_data |>
    arrange(tsne1, tsne2) |>
    mutate(DNAmAge = ifelse(!!sym(clock) > q2, q2, !!sym(clock))) |>
    mutate(DNAmAge = ifelse(DNAmAge < q1, q1, DNAmAge)) |>
    # mutate(clock = ifelse(clock %in% no_accelerat_clocks, clock, paste(clock, "accel"))) |>
    ggplot(aes(x = tsne1, y = tsne2)) +
    geom_point(aes(color = DNAmAge), shape = 19, alpha = 0.6, size = 0.2) +
    labs(color = ifelse(clock %in% no_accelerat_clocks, clock, glue("{clock} accel"))) +
    theme(legend.position = "top") +
    paletteer::scale_color_paletteer_c("pals::coolwarm") +
    better_color_legend

  plotlist[[i]] = plot
}

tot_plot = plotlist |> wrap_plots(nrow = 2) * theme(axis.text.x = element_blank(), axis.text.y = element_blank())
ggsave2("plot/human_dnam_acc_tsne_panels.pdf", tot_plot, w = 40, h = 10)


## Write the plotly obj
library(plotly)
plotly_tsne = function(pp) {
  plot_ly(
    x = pp$tsne1,
    y = pp$tsne2,
    z = pp$tsne3,
    type = "scatter3d",
    hoverinfo = "text",
    text = pp$label,
    marker = list(
      color = pp$DNAmAge,
      colorscale = "Portland",
      alpha = 0.6,
      size = 1,
      showscale = TRUE
    )
  ) |>
    plotly::layout(
      title = clock,
      scene = list(
        title = clock,
        xaxis = list(title = "tsne1"),
        yaxis = list(title = "tsne2"),
        zaxis = list(title = "tsne3")
      ),
      paper_bgcolor = "#000",
      plot_bgcolor = "#000",
      font = list(color = "#fff")
    )
}
annot_df = read_csv("data/human_dnam/annotation.csv")

for (clock in unique(colnames(human_dnam_accer_cl)[2:ncol(human_dnam_accer_cl)])) {
  quants = plot_data[[clock]] |>
    quantile(probs = c(0.05, 0.95), na.rm = T)
  if (quants[[2]] > 10) {
    q1 = round(quants[[1]] / 10) * 10
    q2 = round(quants[[2]] / 10) * 10
  } else {
    q1 = round(quants[[1]] * 2) / 2
    q2 = round(quants[[2]] * 2) / 2
  }

  pp = plot_data |>
    arrange(tsne1, tsne2) |>
    inner_join(annot_df, by = c("sample_id" = "GSM")) |>
    mutate(DNAmAge = ifelse(!!sym(clock) > q2, q2, !!sym(clock))) |>
    mutate(DNAmAge = ifelse(DNAmAge < q1, q1, DNAmAge)) |>
    mutate(label = str_c("</br>Delta Age: ", signif(!!sym(clock), 3), annotation))

  plotly_plot = plotly_tsne(pp)

  saveRDS(plotly_plot, glue("plot/human_dnam_acc_tsne_{clock}.rds"))
  htmlwidgets::saveWidget(plotly_plot, glue("plot/human_dnam_acc_tsne_{clock}.html"), selfcontained = T)
}
