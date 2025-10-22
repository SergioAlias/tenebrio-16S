# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            pcoa.R                                 ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : tenebrio-16S                                     ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-08                                       ║
# ║ Last Modified  : 2025-10-22                                       ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(rlang)
library(magrittr, include.only = "%<>%")
library(tidyverse)
library(qiime2R)
library(patchwork)
library(ggforce)


## Import QIIME 2 files

project_name <- "tenebrio_16S_noC1M1"
local_metadata <- project_name
out <- project_name

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
                    "beta")

jaccard_file_path <- file.path(project_dir,
                               "qiime2/diversity/jaccard_pcoa_results.qza")
bray_curtis_file_path <- file.path(project_dir,
                                   "qiime2/diversity/bray_curtis_pcoa_results.qza")
aitchison_file_path <- file.path(project_dir,
                                 "qiime2/diversity/aitchison_pcoa_results.qza")

jaccard <- read_qza(jaccard_file_path)
bray_curtis <- read_qza(bray_curtis_file_path)
aitchison <- read_qza(aitchison_file_path)

metadata <- read.csv(file.path("/home/sergio/scratch",
                               local_metadata,
                               "metadata.tsv"),
                     sep = "\t")

colnames(metadata)[1] <- "SampleID"

## Get variance explained by each PCo 

jaccard_pco1 <- round(jaccard[["data"]]$ProportionExplained$PC1 * 100, 2)
jaccard_pco2 <- round(jaccard[["data"]]$ProportionExplained$PC2 * 100, 2)

bray_curtis_pco1 <- round(bray_curtis[["data"]]$ProportionExplained$PC1 * 100, 2)
bray_curtis_pco2 <- round(bray_curtis[["data"]]$ProportionExplained$PC2 * 100, 2)

aitchison_pco1 <- round(aitchison[["data"]]$ProportionExplained$PC1 * 100, 2)
aitchison_pco2 <- round(aitchison[["data"]]$ProportionExplained$PC2 * 100, 2)


## Colors and shapes

source("/home/sergio/projects/tenebrio_16S_noC1M1/colors.R")

## PCoA plots

my_breaks <- c("Control", "DON", "AFB1", "FB1") 
my_plotmath_labels <- c("Control", "DON", bquote(""*AFB[1]), bquote(""*FB[1]))

### Jaccard (confidence ellipses)

p_j <- jaccard$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Mycotoxin`, shape = `Mycotoxin`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = beta_colors, name = "Mycotoxin") +
  scale_shape_manual(values = beta_shapes, name = "Mycotoxin") +
  xlab(paste0("PCo-1 | ", jaccard_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", jaccard_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Mycotoxin`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

pdf(file.path(outdir, "pcoa_jaccard.pdf"),
    height = 5,
    width = 7)

p_j

dev.off()

### Bray-Curtis (hull and centroids)

plot_data_bray <- bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata)

centroids_bray <- plot_data_bray %>%
  group_by(`Mycotoxin`) %>%
  summarise(
    PC1_cent = mean(PC1),
    PC2_cent = mean(PC2),
    .groups = 'drop'
  )

hull_data_bray <- plot_data_bray %>%
  group_by(`Mycotoxin`) %>%
  slice(chull(PC1, PC2)) 

plot_data_with_centroids <- plot_data_bray %>%
  left_join(centroids_bray, by = "Mycotoxin")

p_b <- ggplot(plot_data_with_centroids,
              aes(x = PC1, y = PC2, color = `Mycotoxin`, shape = `Mycotoxin`)) +
  geom_polygon(data = hull_data_bray, 
               aes(fill = `Mycotoxin`), 
               alpha = 0.2, 
               show.legend = FALSE,
               linetype = "dotted",
               linewidth = 0.5) +
  geom_segment(aes(xend = PC1_cent, yend = PC2_cent),
               alpha = 0.5, 
               show.legend = FALSE) +
  geom_point(alpha = 1, size = 3) + 
  geom_point(data = centroids_bray, 
             aes(x = PC1_cent, y = PC2_cent),
             size = 2, 
             shape = 5,
             show.legend = FALSE) +
  theme_bw() +
  scale_color_manual(values = beta_colors, name = NULL,
                     breaks = my_breaks,
                     labels = my_plotmath_labels) +
  scale_shape_manual(values = beta_shapes, name = NULL,
                     breaks = my_breaks,
                     labels = my_plotmath_labels) +
scale_fill_manual(values = beta_colors, name = "Mycotoxin") + 
  xlab(paste0("PCo-1 | ", bray_curtis_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", bray_curtis_pco2, "% of variance explained")) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

pdf(file.path(outdir, "pcoa_bray_curtis.pdf"),
    height = 5,
    width = 7)

p_b

dev.off()

### Aitchison (confidence ellipses)

p_a <- aitchison$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Mycotoxin`, shape = `Mycotoxin`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = beta_colors, name = "Mycotoxin") +
  scale_shape_manual(values = beta_shapes, name = "Mycotoxin") +
  xlab(paste0("PCo-1 | ", aitchison_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", aitchison_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Mycotoxin`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

pdf(file.path(outdir, "pcoa_aitchison.pdf"),
    height = 5,
    width = 7)

p_a

dev.off()
