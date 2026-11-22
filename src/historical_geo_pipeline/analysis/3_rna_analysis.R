library(tidyverse)
library(glue)
library(patchwork)
library(recipes)
library(embed)
library(plotly)

# Load data
human_rna_age <- read_csv("result/human_rna_age.csv")
human_rna_meta <- read_csv("result/gpt4_human_rna_metadata.csv")
mouse_rna_age <- read_csv("result/mouse_rna_age.csv")
mouse_rna_meta <- read_csv("result/gpt4_mouse_rna_metadata.csv")

human_umap = human_rna_age %>%
  recipes::recipe(formula(~.), data = select(human_rna_age, gse_id, gsm_id, everything(), -ends_with("_std"))) %>%
  recipes::step_impute_mean(recipes::all_numeric_predictors()) %>%
  recipes::step_normalize(recipes::all_numeric_predictors()) %>%
  embed::step_umap(recipes::all_numeric_predictors(), neighbors = 200, num_comp = 3) %>%
  recipes::prep() %>%
  recipes::juice()

write_csv(human_umap, "result/human_rna_umap_age.csv")

mouse_umap = mouse_rna_age %>%
  recipes::recipe(formula(~.), data = select(mouse_rna_age, gse_id, gsm_id, everything(), -ends_with("_std"))) %>%
  recipes::step_impute_mean(recipes::all_numeric_predictors()) %>%
  recipes::step_normalize(recipes::all_numeric_predictors()) %>%
  embed::step_umap(recipes::all_numeric_predictors(), neighbors = 200, num_comp = 3) %>%
  recipes::prep() %>%
  recipes::juice()

write_csv(mouse_umap, "result/mouse_rna_umap_age.csv")

# Function to process and analyze data
analyze_rna_data <- function(age_data, meta_data, species) {
  # Merge age and metadata
  data <- age_data %>%
    inner_join(meta_data, by = c("gse_id", "gsm_id"))
  
  # Convert age to numeric and filter valid ages
  data <- data %>%
    mutate(age = as.numeric(age)) %>%
    filter(age >= 0 & age <= ifelse(species == "human", 120, 3))
  
  # Get clock names
  clocks <- setdiff(names(data), c("gse_id", "gsm_id", "age", names(meta_data)))
  
  # Perform UMAP
  set.seed(123)
  umap_recipe <- recipe(~., data = data %>% select(all_of(clocks))) %>%
    step_impute_knn(all_predictors()) %>%
    step_normalize(all_predictors()) %>%
    step_umap(all_predictors(), num_comp = 3, neighbors = 15)
  
  umap_prep <- prep(umap_recipe)
  umap_data <- juice(umap_prep)
  
  umap_data <- umap_data %>%
    bind_cols(data %>% select(gse_id, gsm_id, age))
  
  # Plot UMAP for each clock
  plotlist <- list()
  
  for (clock in clocks) {
    plot_data <- umap_data %>%
      mutate(clock_value = data[[clock]])
    
    p <- plot_data %>%
      ggplot(aes(x = UMAP1, y = UMAP2, color = clock_value)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis_c() +
      labs(title = glue("{species} RNA - {clock}"), color = clock) +
      theme_minimal()
    
    plotlist[[clock]] <- p
    
    ggsave(glue("plot/{species}_rna_{clock}_umap.pdf"), p, width = 10, height = 8)
    
    # Create interactive 3D plot
    plot_3d <- plot_ly(
      data = plot_data,
      x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
      color = ~clock_value,
      type = "scatter3d", mode = "markers",
      marker = list(size = 3, opacity = 0.8)
    ) %>%
      layout(
        scene = list(
          xaxis = list(title = "UMAP1"),
          yaxis = list(title = "UMAP2"),
          zaxis = list(title = "UMAP3")
        ),
        title = glue("{species} RNA - {clock}")
      )
    
    saveRDS(plot_3d, glue("plot/{species}_rna_{clock}_umap_3d.rds"))
    htmlwidgets::saveWidget(plot_3d, glue("plot/{species}_rna_{clock}_umap_3d.html"), selfcontained = TRUE)
  }
  
  # Combine all UMAP plots
  combined_plot <- wrap_plots(plotlist, ncol = 3)
  ggsave(glue("plot/{species}_rna_all_clocks_umap.pdf"), combined_plot, width = 20, height = 15)
  
  return(list(umap_data = umap_data, plotlist = plotlist))
}

# Analyze human RNA data
human_results <- analyze_rna_data(human_rna_age, human_rna_meta, "human")

# Analyze mouse RNA data
mouse_results <- analyze_rna_data(mouse_rna_age, mouse_rna_meta, "mouse")
