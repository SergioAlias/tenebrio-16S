# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            rarefaction.R                          ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : tenebrio-16S                                     ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2026-06-04                                       ║
# ║ Last Modified  : 2026-06-08                                       ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(tidyverse)
library(magrittr, include.only = "%<>%")
library(patchwork)

## Colors and shapes

source("/home/salias/projects/tenebrio-16S/colors.R")

## Functions

plot_rarefaction <- function(csv_file, metadata, y_title, threshold = NULL) {
  data <- read_csv(csv_file, show_col_types = FALSE)
  
  long_data <- data %>%
    pivot_longer(
      cols = -`sample-id`, 
      names_to = c("depth", "iteration"),
      names_pattern = "depth-(.*)_iter-(.*)",
      values_to = "observed_features",
      values_drop_na = TRUE
    ) %>%
    mutate(
      depth = as.numeric(depth),
      iteration = as.numeric(iteration)
    )
  
  summary_data <- long_data %>%
    group_by(`sample-id`, depth) %>%
    summarise(
      mean_features = mean(observed_features),
      .groups = 'drop'
    )
  
  summary_data <- summary_data %>%
    left_join(metadata, by = c("sample-id" = "SampleID"))
  
  legend_order <- c("Control", "DON", "AFB1", "FB1") 
  legend_labels <- list("Control", "DON", bquote(AFB[1]), bquote(FB[1]))
  
  p <- ggplot(summary_data, aes(x = depth, y = mean_features, group = `sample-id`, color = Mycotoxin)) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    scale_color_manual(
      values = beta_colors,
      breaks = legend_order,
      labels = legend_labels
    ) + 
    labs(
      x = "Sequencing depth",
      y = y_title,            
      color = "Mycotoxin"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(t = 20, r = 5, b = 5, l = 5) 
    )
  
  if (!is.null(threshold)) {
    p <- p + geom_vline(
      xintercept = threshold, 
      linetype = "dotted", 
      color = "black", 
      linewidth = 0.8,
      alpha = 0.8
    ) +
      annotate(
        "text", 
        x = threshold, 
        y = Inf,                         
        label = as.character(threshold), 
        angle = 0,                        
        vjust = -1,                      
        hjust = 0.5,
        size = 3.5                        
      ) +
      coord_cartesian(clip = "off")       
  }
  return(p)
}

clean_metadata <- function(path) {
  df <- read.csv(path, sep = "\t")
  df %<>%
    rename(SampleID = 1)
  return(df)
}

## Import files

out <- "tenebrio-16S"
outdir <- file.path("/home/salias/scratch",
                    out,
                    "raref")
metadata <- clean_metadata(file.path("/home/salias/scratch",
                                      out,
                                      "metadata.tsv"))



## Rarefaction plots

of_16s <- plot_rarefaction(file.path("/home/salias/scratch", 
                                     out, 
                                     "raref/observed_features.csv"),
                           metadata,
                           y_title = "Observed features",
                           threshold = 64351)

s_16s <- plot_rarefaction(file.path("/home/salias/scratch", 
                                    out, 
                                    "raref/shannon.csv"),
                          metadata,
                          y_title = "Shannon",
                          threshold = 64351)


p_raref <- of_16s / plot_spacer() / s_16s + 
  plot_layout(
    guides = "collect",
    heights = c(1, 0.1, 1) 
  ) & 
  theme(
    legend.position = "right"
  )


## Save plots

pdf(file.path(outdir, "rarefaction_curves.pdf"),
    width = 9,
    height = 7)

p_raref

dev.off()

png(file.path(outdir, "rarefaction_curves.png"),
    width = 9,
    height = 7,
    units = "in",
    res = 300)

p_raref

dev.off()
