# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            taxonomy.R                             ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : tenebrio-16S                                     ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-06                                       ║
# ║ Last Modified  : 2025-10-28                                       ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(tidyverse)
library(ggplot2)
library(qiime2R)


## Import QIIME 2 files

project_name <- "tenebrio_16S_noC1M1"
local_metadata <- "tenebrio_16S_noC1M1"
out <- "tenebrio_16S_noC1M1"
source("/home/sergio/projects/tenebrio_16S_noC1M1/colors.R")

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")
project_dir <- file.path(cluster_path,
                         "scratch/salias/projects",
                         project_name)
outdir <- file.path("/home/sergio/scratch",
                    out,
                    "taxonomy")


dada2_file_path <- file.path(project_dir,
                             "qiime2/feature_tables/filtered_table.qza")
metadata_file_path <- file.path("/home/sergio/scratch",
                                local_metadata,
                                "metadata.tsv")
taxonomy_file_path <- file.path(project_dir,
                                "qiime2/taxonomy/taxonomy.qza")


## Custom barplot

df <- read_qza(dada2_file_path)$data
df_long <- data.frame(
  row = rep(rownames(df), times = ncol(df)),
  column = rep(colnames(df), each = nrow(df)),
  value = as.vector(df)
)

tax <- read_qza(taxonomy_file_path)$data %>% parse_taxonomy()

tax$row <- rownames(tax)
rownames(tax) <- NULL

df_long %<>%
  left_join(tax %>% select(row, Genus), by = "row")


metadata <- read.csv(metadata_file_path, sep = "\t", header = TRUE)
metadata$column <- metadata$ID
metadata$ID <- NULL

df_long %<>%
  left_join(metadata %>% select(column, Mycotoxin), by = "column")

df_long %<>%
  mutate(Genus = replace_na(Genus, "Unassigned"))

top_14_genera <- df_long %>%
  filter(Genus != "Unassigned" & !grepl("Incertae_sedis", Genus)) %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(value, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 14) %>%
  pull(Genus)

df_long %<>%
  mutate(
    Genus_plot = case_when(
      Genus %in% top_14_genera ~ Genus,
      Genus == "Unassigned" ~ "Unassigned",
      TRUE ~ "Other"
    )
  )

plot_levels <- c(top_14_genera, "Other", "Unassigned")
df_long$Genus_plot <- factor(df_long$Genus_plot, levels = plot_levels)

plot_colors <- setNames(
  c(barplot_colors[1:15], "grey80"),
  plot_levels
)

legend_labels <- sapply(plot_levels, function(label) {
  if (label %in% c("Unassigned", "Other")) {
    return(label) # Sin cursiva
  } else {
    return(bquote(italic(.(label))))
  }
})

create_custom_barplot <- function(vertical = FALSE) {
  
  myco_order <- c("Control", "DON", "AFB1", "FB1")
  myco_labels <- c("Control", "DON", "AFB[1]", "FB[1]")
  
  df_long$Mycotoxin <- factor(
    df_long$Mycotoxin,
    levels = myco_order,
    labels = myco_labels
  )
  
  p <- ggplot(df_long) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(
      name = "Genus",
      values = plot_colors,
      breaks = plot_levels,
      labels = legend_labels
    ) +
    theme(
      text = element_text(size = 15),
      legend.title = element_text(size = 15),
      legend.position = "right",
      panel.background = element_blank(),
      panel.border = element_blank(),
      strip.text = element_text(
        size = 15,
        color = "black",
        hjust = 0.5,
        vjust = 0.5
      ),
      strip.placement = "outside",
      strip.background = element_rect(fill = NA, color = NA)
    ) +
    guides(fill = guide_legend(nrow = 16, ncol = 1))
  
  if (vertical) {
    p <- p +
      aes(fill = Genus_plot, y = value, x = column) +
      labs(y = "Relative abundance", x = NULL) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12, color = "black", margin = margin(r = -2)),
        axis.ticks.y = element_blank()
      ) +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
      scale_x_discrete(expand = c(0, 0.7)) +
      facet_grid(. ~ Mycotoxin, switch = "x", scale = "free_x", labeller = label_parsed)
  } else {
    p <- p +
      aes(fill = Genus_plot, x = value, y = column) +
      labs(x = "Relative abundance", y = NULL) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 12, color = "black", margin = margin(r = -2)),
        axis.ticks.x = element_blank()
      ) +
      scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
      scale_y_discrete(expand = c(0, 0.7)) +
      facet_grid(Mycotoxin ~ ., switch = "y", scale = "free_y", labeller = label_parsed)
  }
  
  return(p)
}

hor <- create_custom_barplot()
ver <- create_custom_barplot(vertical = TRUE)

## Save plots

pdf(file.path(outdir, "custom_barplot_horizontal.pdf"),
    width = 9,
    height = 7)

hor

dev.off()

png(file.path(outdir, "custom_barplot_horizontal.png"),
    width = 9,
    height = 7,
    units = "in",
    res = 300)

hor

dev.off()

pdf(file.path(outdir, "custom_barplot_vertical.pdf"),
    width = 9,
    height = 7)

ver

dev.off()

png(file.path(outdir, "custom_barplot_vertical.png"),
    width = 9,
    height = 7,
    units = "in",
    res = 300)

ver

dev.off()
