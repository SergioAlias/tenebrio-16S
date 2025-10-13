# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            alpha.R                                ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : tenebrio-16S                                     ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-08                                       ║
# ║ Last Modified  : 2025-10-08                                       ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(tidyverse)
library(qiime2R)
library(ggpubr)
library(patchwork)
library(rstatix)
library(multcompView)
library(patchwork)


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
                    "alpha")

metadata <- read.csv(file.path("/home/sergio/scratch",
                               local_metadata,
                               "metadata.tsv"),
                     sep = "\t")

colnames(metadata)[1] <- "SampleID"

shannon_file_path <- file.path(project_dir,
                               "qiime2/diversity/shannon_vector.qza")
simpson_file_path <- file.path(project_dir,
                               "qiime2/diversity/simpson_vector.qza")
chao1_file_path <- file.path(project_dir,
                             "qiime2/diversity/chao1_vector.qza")

shannon <- read_qza(shannon_file_path)
shannon <- shannon$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(shannon)

simpson <- read_qza(simpson_file_path)
simpson <- simpson$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(simpson)

chao1 <- read_qza(chao1_file_path)
chao1 <- chao1$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(chao1)
metadata$Mycotoxin <- factor(metadata$Mycotoxin, levels = c("Control", "DON", "AFB1", "FB1"))

## Colors and shapes

source("/home/sergio/projects/tenebrio_16S_noC1M1/colors.R")

## Comparisons

comparisons <- combn(unique(metadata$Mycotoxin), 2, simplify = FALSE)

comparisons <- lapply(setdiff(unique(metadata$Mycotoxin), "Control"), function(x) c("Control", x))


## Alpha boxplots

shannon_pos_stat <- 3.7
shannon_pos_stat_paired <- c(3.4, 3.3, 3.5)
simpson_pos_stat <- 0.91
simpson_pos_stat_paired <- c(0.84, 0.82, 0.86)
chao1_pos_stat <- 115
chao1_pos_stat_paired <- c(106.5, 104, 109)

### Shannon

shannon <- metadata %>%
  ggboxplot("Mycotoxin", "shannon_entropy",
            color = "Mycotoxin",
            fill = "Mycotoxin",
            alpha = 0.1,
            palette = alpha_colors,
            add = "jitter",
            shape = "Mycotoxin") +
  theme(legend.position = "none") +
  scale_shape_manual(values = alpha_shapes, name = "Mycotoxin") +
  ylab("Shannon") +
  stat_compare_means(label.y = shannon_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons,
                     label.y = shannon_pos_stat_paired)

pdf(file.path(outdir, "shannon.pdf"))

shannon

dev.off()

### Inverse Simpson

simpson <- metadata %>%
  ggboxplot("Mycotoxin", "simpson",
            color = "Mycotoxin",
            fill = "Mycotoxin",
            alpha = 0.1,
            palette = alpha_colors,
            add = "jitter",
            shape = "Mycotoxin") +
  theme(legend.position = "none") +
  scale_shape_manual(values = alpha_shapes, name = "Mycotoxin") +
  ylab("Inverse Simpson") +
  stat_compare_means(label.y = simpson_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons,
                     label.y = simpson_pos_stat_paired)

pdf(file.path(outdir, "simpson.pdf"))

simpson

dev.off()

### Chao1

chao1 <- metadata %>%
  ggboxplot("Mycotoxin", "chao1",
            color = "Mycotoxin",
            fill = "Mycotoxin",
            alpha = 0.1,
            palette = alpha_colors,
            add = "jitter",
            shape = "Mycotoxin") +
  theme(legend.position = "none") +
  scale_shape_manual(values = alpha_shapes, name = "Mycotoxin") +
  ylab("Chao1") +
  stat_compare_means(label.y = chao1_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons,
                     label.y = chao1_pos_stat_paired)

pdf(file.path(outdir, "chao1.pdf"))

chao1

dev.off()

### Grouped plot

pdf(file.path(outdir, "alpha_patched.pdf"),
    width = 10,
    height = 6)

(chao1 + theme(legend.position="none") +
    shannon + theme(legend.position="none") +
    simpson + theme(legend.position="none") +
    #theme(plot.tag.position = "topleft")) +
  plot_layout(axis_titles = "collect") )#+
  #plot_annotation(tag_levels = 'A')

dev.off()
