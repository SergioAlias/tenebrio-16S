# ╔═══════════════════════════════════════════════════════════════════╗
# ║                       abundance_tables.R                          ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : tenebrio-16S                                     ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-07                                       ║
# ║ Last Modified  : 2025-10-07                                       ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(qiime2R)
library(tidyverse)
library(openxlsx)

## Table construction

project_name <- "tenebrio_16S_noC1M1"
local_metadata <- "tenebrio_16S_noC1M1"
out <- "tenebrio_16S_noC1M1"

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
                    "feature_table")


dada2_file_path <- file.path(project_dir,
                             "qiime2/feature_tables/filtered_table.qza")
metadata_file_path <- file.path("/home/sergio/scratch",
                                local_metadata,
                                "metadata.tsv")
taxonomy_file_path <- file.path(project_dir,
                                "qiime2/taxonomy/taxonomy.qza")

df <- read_qza(dada2_file_path)$data
tax <- read_qza(taxonomy_file_path)$data %>% parse_taxonomy()
merged_df <- merge(df, tax, by = "row.names")

write.xlsx(merged_df, file = file.path(outdir, "feature_table.xlsx"))
