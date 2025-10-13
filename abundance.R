# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            abundance.R                            ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : tenebrio-16S                                     ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-08                                       ║
# ║ Last Modified  : 2025-10-13                                       ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

# Libraries

library(magrittr, include.only = "%<>%")
library(qiime2R)
library(readr)
library(tidyverse)
library(EnhancedVolcano)
library(patchwork)

## Colors and shapes

source("/home/sergio/projects/tenebrio_16S_noC1M1/colors.R")

## Functions

source("/home/sergio/projects/diversity-cereal/parse_ancombc.R")

volcanoFromAncombc <- function(qza_path,
                               log2fc_col,
                               pval_col,
                               up_color,
                               down_color,
                               up_shape,
                               down_shape,
                               up_legend,
                               down_legend,
                               ...,
                               lab_col = NA,
                               log2fc_cutoff = 2,
                               pval_cutoff = 0.05,
                               ns_color = "grey45",
                               ns_shape = 4,
                               ns_legend = "NS",
                               taxonomy_df = taxonomy)
{
  ancombc <- import_ancombc(qza_path)
  ancombc %<>% left_join(taxonomy_df)
  ancombc %<>% mutate(
    labnames = coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom),
    labnames = if_else(
      labnames %in% c(Species, Genus),
      paste0("italic('", labnames, "')"),
      paste0("'", labnames, "'")),
    labnames = gsub("_", " ", labnames))
  
  keyvals_col <- ifelse(
    ancombc[[log2fc_col]] < -log2fc_cutoff & ancombc[[pval_col]] < pval_cutoff, down_color,
    ifelse(ancombc[[log2fc_col]] > log2fc_cutoff & ancombc[[pval_col]] < pval_cutoff, up_color,
           ns_color))
  keyvals_col[is.na(keyvals_col)] <- ns_color
  names(keyvals_col)[keyvals_col == up_color] <- up_legend
  names(keyvals_col)[keyvals_col == ns_color] <- ns_legend
  names(keyvals_col)[keyvals_col == down_color] <- down_legend
  
  keyvals_shape <- keyvals_col
  keyvals_shape[keyvals_shape == up_color] <- up_shape
  keyvals_shape[keyvals_shape == ns_color] <- ns_shape
  keyvals_shape[keyvals_shape == down_color] <- down_shape
  keyvals_shape %<>% as.integer()
  names(keyvals_shape) <- names(keyvals_col)
  
  lab_arg <- if (is.na(lab_col)) NA else ancombc[[lab_col]]
  
  v_plot <- ancombc %>%
    EnhancedVolcano(lab = lab_arg,
                    x = {{log2fc_col}},
                    y = {{pval_col}},
                    pCutoff = pval_cutoff,
                    FCcutoff = log2fc_cutoff,
                    colCustom = keyvals_col,
                    shapeCustom = keyvals_shape,
                    ...) +
    guides(color = guide_legend("Combined Legend",
                                override.aes = list(alpha=1)),
           shape = guide_legend("Combined Legend")) +
    theme_classic() +
    theme(legend.title=element_blank(),
          legend.position="top")
  
  return(v_plot)
}

## Import QIIME 2 files

project_name <- "tenebrio_16S_noC1M1"
out <- project_name
afl_tag <- "AFB1"
con_tag <- "Control"
don_tag <- "DON"
fum_tag <- "FB1"

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
                    "abundance")

treatment_con_file_path <- file.path(project_dir,
                                     paste0("qiime2/abundance/Mycotoxin_",
                                            con_tag,
                                            "/filtered_ancombc.qza"))

taxonomy_file_path <- file.path(project_dir,
                                "qiime2/taxonomy/taxonomy.qza")

taxonomy <- read_qza(taxonomy_file_path)$data

taxonomy %<>% parse_taxonomy() %>% rownames_to_column("id")


## Volcano plots

### AFB1 vs CON

v_AFB1_vs_CON <- volcanoFromAncombc(qza_path = treatment_con_file_path,
                                    lab_col = "labnames",
                                    log2fc_col = paste0("Mycotoxin", afl_tag, "_lfc"),
                                    pval_col = paste0("Mycotoxin", afl_tag, "_q_val"),
                                    up_color = da_colors[["AFB1"]],
                                    down_color = da_colors[["Control"]],
                                    up_shape = da_shapes[["AFB1"]],
                                    down_shape = da_shapes[["Control"]],
                                    up_legend = "DA (AFB1)",
                                    down_legend = "DA (CON)",
                                    colAlpha = 1,
                                    ylab = bquote(~Log[10]~ "Q-value"),
                                    title = NULL,
                                    subtitle = NULL,
                                    caption = NULL,
                                    xlim = c(-5, 5),
                                    ylim = c(0, 8),
                                    drawConnectors = TRUE,
                                    typeConnectors = "closed",
                                    widthConnectors = 0,
                                    pointSize = 4.0,
                                    labSize = 3.5,
                                    boxedLabels = TRUE,
                                    parseLabels = TRUE)

pdf(file.path(outdir, "volcano_AFB1_vs_CON.pdf"))

v_AFB1_vs_CON

dev.off()

### DON vs CON

v_DON_vs_CON <- volcanoFromAncombc(qza_path = treatment_con_file_path,
                                    lab_col = "labnames",
                                    log2fc_col = paste0("Mycotoxin", don_tag, "_lfc"),
                                    pval_col = paste0("Mycotoxin", don_tag, "_q_val"),
                                    up_color = da_colors[["DON"]],
                                    down_color = da_colors[["Control"]],
                                    up_shape = da_shapes[["DON"]],
                                    down_shape = da_shapes[["Control"]],
                                    up_legend = "DA (DON)",
                                    down_legend = "DA (CON)",
                                    colAlpha = 1,
                                    ylab = bquote(~Log[10]~ "Q-value"),
                                    title = NULL,
                                    subtitle = NULL,
                                    caption = NULL,
                                    xlim = c(-5, 5),
                                    ylim = c(0, 8),
                                    drawConnectors = TRUE,
                                    typeConnectors = "closed",
                                    widthConnectors = 0,
                                    pointSize = 4.0,
                                    labSize = 3.5,
                                    boxedLabels = TRUE,
                                    parseLabels = TRUE)

pdf(file.path(outdir, "volcano_DON_vs_CON.pdf"))

v_DON_vs_CON

dev.off()

### FB1 vs CON

v_FB1_vs_CON <- volcanoFromAncombc(qza_path = treatment_con_file_path,
                                   lab_col = "labnames",
                                   log2fc_col = paste0("Mycotoxin", fum_tag, "_lfc"),
                                   pval_col = paste0("Mycotoxin", fum_tag, "_q_val"),
                                   up_color = da_colors[["FB1"]],
                                   down_color = da_colors[["Control"]],
                                   up_shape = da_shapes[["FB1"]],
                                   down_shape = da_shapes[["Control"]],
                                   up_legend = "DA (FB1)",
                                   down_legend = "DA (CON)",
                                   colAlpha = 1,
                                   ylab = bquote(~Log[10]~ "Q-value"),
                                   title = NULL,
                                   subtitle = NULL,
                                   caption = NULL,
                                   xlim = c(-5, 5),
                                   ylim = c(0, 8),
                                   drawConnectors = TRUE,
                                   typeConnectors = "closed",
                                   widthConnectors = 0,
                                   pointSize = 4.0,
                                   labSize = 3.5,
                                   boxedLabels = TRUE,
                                   parseLabels = TRUE)

pdf(file.path(outdir, "volcano_FB1_vs_CON.pdf"))

v_FB1_vs_CON

dev.off()

### Grouped plots

pdf(file.path(outdir, "patched_volcano.pdf"),
    width = 15,
    height = 8)

(v_AFB1_vs_CON +
    v_DON_vs_CON +
    v_FB1_vs_CON +
    theme(plot.tag.position = "topleft")) +
  plot_layout(axis_titles = "collect") # + # ,
# guides = "collect") +
# plot_annotation(tag_levels = 'A')

dev.off()
