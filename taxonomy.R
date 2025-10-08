# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            taxonomy.R                             ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : tenebrio-16S                                     ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-06                                       ║
# ║ Last Modified  : 2025-10-08                                       ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(tidyverse)
library(file2meco)
library(microeco)
library(ggplot2)
library(ggnested)
library(ggh4x)
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

meco <- qiime2meco(dada2_file_path,
                   sample_table = metadata_file_path,
                   taxonomy_table = taxonomy_file_path)


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
  group_by(Genus) %>%
  mutate(total_abundance = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Genus = reorder(Genus, -total_abundance))

top <- c(levels(df_long$Genus)[1:15], NA)

legend_labels <- sapply(top, function(label) {
  if (is.na(label)) {
    "Unassigned"
  } else {
    bquote(italic(.(label)))
  }
})

custom_barplot <- ggplot(df_long, aes(fill=Genus, y=value, x=column)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(
    values = setNames(c(barplot_colors, "grey"), top),
    breaks = top,
    labels = legend_labels
  ) +
  labs(x = "", y = "Relative abundance (%)") +
  theme(text = element_text(size = 15),
        legend.title = element_text(size = 15),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 14, margin = margin(r =-2)),
        axis.ticks.y = element_blank(), 
        strip.text = element_text(
          size = 15,                       
          face = "bold",                   
          color = "black",                 
          hjust = 0.5,                     
          vjust = 0.5                      
        ),
        strip.background = element_rect(
          fill = NA,
          color = NA
        )) + 
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) x * 100) +
  scale_x_discrete(expand = c(0, 0.7)) +
  facet_grid(~Mycotoxin, scale = "free_x", space = "free_x")

pdf(file.path(outdir, "custom_barplot.pdf"),
    width = 9)

custom_barplot

dev.off()

## Relabel UNITE prefixes for cleaner plotting

meco$tax_table$Phylum <- gsub("p__", "", meco$tax_table$Phylum)
meco$tax_table$Family <- gsub("f__", "", meco$tax_table$Family)

# Create trans_abund objects and plot stuff

## Nested barplot (Phylum / Class)

t_stacked_phylum <- trans_abund$new(dataset = meco,
                                    taxrank = "Class",
                                    ntaxa = 20,
                                    delete_taxonomy_prefix = TRUE,
                                    high_level = "Phylum",
                                    prefix = "c__")

pdf(file.path(outdir, "barplot_class.pdf"),
    width = 9)

t_stacked_phylum$plot_bar(ggnested = TRUE,
                          high_level_add_other = TRUE,
                          xtext_keep = FALSE,
                          # xtext_angle = 90,
                          # xtext_size = 6,
                          facet = c("Mycotoxin"),
                          others_color = "grey90") + 
  theme(ggh4x.facet.nestline = element_line(colour = "grey95"))

dev.off()

## Nested barplot (Family / Genus)

t_stacked_family <- trans_abund$new(dataset = meco,
                                    taxrank = "Genus",
                                    ntaxa = 15,
                                    delete_taxonomy_prefix = TRUE,
                                    high_level = "Family",
                                    prefix = "g__")

pdf(file.path(outdir, paste0(color_palette, "_barplot_genus.pdf")),
    width = 7.5,
    height = 5.5)

t_stacked_family$plot_bar(ggnested = TRUE,
                          high_level_add_other = TRUE,
                          xtext_keep = FALSE,
                          color_values = barplot_colors,
                          # xtext_angle = 90,
                          # xtext_size = 6,
                          facet = c("Mycotoxin")) + 
  theme(ggh4x.facet.nestline = element_line(colour = "grey95"))

dev.off()

