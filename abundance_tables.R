# ╔═══════════════════════════════════════════════════════════════════╗
# ║                       abundance_tables.R                          ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : tenebrio-16S                                     ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-09-26                                       ║
# ║ Last Modified  : 2025-09-26                                       ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(qiime2R)
library(tidyverse)
library(openxlsx)

## Functions

source("/home/sergio/projects/diversity-cereal/parse_ancombc.R")

## Import QIIME 2 files

project_name <- "tenebrio_16S_noC1M1"
out <- "tenebrio_16S_noC1M1"
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

myco_con_file_path <- file.path(project_dir,
                                     paste0("qiime2/abundance/Mycotoxin_",
                                            con_tag,
                                            "/filtered_ancombc.qza"))
myco_afl_file_path <- file.path(project_dir,
                                     paste0("qiime2/abundance/Mycotoxin_",
                                            afl_tag,
                                            "/filtered_ancombc.qza"))
myco_don_file_path <- file.path(project_dir,
                                paste0("qiime2/abundance/Mycotoxin_",
                                       don_tag,
                                       "/filtered_ancombc.qza"))

isc_con_file_path <- file.path(project_dir,
                                paste0("qiime2/abundance/IsControl_",
                                       con_tag,
                                       "/filtered_ancombc.qza"))

myco_con <- import_ancombc(myco_con_file_path)
myco_afl <- import_ancombc(myco_afl_file_path)
myco_don <- import_ancombc(myco_don_file_path)
isc_con <- import_ancombc(isc_con_file_path)

taxonomy_file_path <- file.path(project_dir,
                                "qiime2/taxonomy/taxonomy.qza")

taxonomy <- read_qza(taxonomy_file_path)$data

taxonomy_parsed <- taxonomy %>% parse_taxonomy() %>% rownames_to_column("id")

colnames(taxonomy)[1] <- "id"

taxonomy <- taxonomy[, c("id", "Confidence")]

taxonomy_parsed %<>% left_join(taxonomy)
myco_con %<>% left_join(taxonomy_parsed)
myco_afl %<>% left_join(taxonomy_parsed)
myco_don %<>% left_join(taxonomy_parsed)
isc_con %<>% left_join(taxonomy_parsed)

## Table filtering

qval_thr <- 0.05
lfc_thr <- 2

fum_tag <- "MycotoxinFB1"

AFB1vsCON <- myco_con[, !grepl(don_tag, colnames(myco_con)) & !grepl(fum_tag, colnames(myco_con))]
DONvsCON <- myco_con[, !grepl(afl_tag, colnames(myco_con)) & !grepl(fum_tag, colnames(myco_con))]
FB1vsCON <- myco_con[, !grepl(afl_tag, colnames(myco_con)) & !grepl(don_tag, colnames(myco_con))]
DONvsAFB1 <- myco_afl[, !grepl(con_tag, colnames(myco_afl)) & !grepl(fum_tag, colnames(myco_afl))]
FB1vsAFB1 <- myco_afl[, !grepl(con_tag, colnames(myco_afl)) & !grepl(don_tag, colnames(myco_afl))]
FB1vsDON <- myco_don[, !grepl(con_tag, colnames(myco_don)) & !grepl(afl_tag, colnames(myco_don))]
MYCOvsCON <- isc_con

fum_tag <- "FB1"

AFB1vsCON_up <- AFB1vsCON %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", afl_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc >= lfc_thr) %>%
  arrange(desc(lfc))

AFB1vsCON_down <- AFB1vsCON %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", afl_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc <= -lfc_thr) %>%
  arrange(lfc)

DONvsCON_up <- DONvsCON %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", don_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc >= lfc_thr) %>%
  arrange(desc(lfc))

DONvsCON_down <- DONvsCON %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", don_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc <= -lfc_thr) %>%
  arrange(lfc)

FB1vsCON_up <- FB1vsCON %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", fum_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc >= lfc_thr) %>%
  arrange(desc(lfc))

FB1vsCON_down <- FB1vsCON %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", fum_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc <= -lfc_thr) %>%
  arrange(lfc)

DONvsAFB1_up <- DONvsAFB1 %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", don_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc >= lfc_thr) %>%
  arrange(desc(lfc))

DONvsAFB1_down <- DONvsAFB1 %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", don_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc <= -lfc_thr) %>%
  arrange(lfc)

FB1vsAFB1_up <- FB1vsAFB1 %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", fum_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc >= lfc_thr) %>%
  arrange(desc(lfc))

FB1vsAFB1_down <- FB1vsAFB1 %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", fum_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc <= -lfc_thr) %>%
  arrange(lfc)

FB1vsDON_up <- FB1vsDON %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", fum_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc >= lfc_thr) %>%
  arrange(desc(lfc))

FB1vsDON_down <- FB1vsDON %>%
  rename_with(~str_remove(.x, paste0("Mycotoxin", fum_tag, "_"))) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc <= -lfc_thr) %>%
  arrange(lfc)

MYCOvsCON_up <- MYCOvsCON %>%
  rename_with(~str_remove(.x, "IsControlMycotoxin_")) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc >= lfc_thr) %>%
  arrange(desc(lfc))

MYCOvsCON_down <- MYCOvsCON %>%
  rename_with(~str_remove(.x, "IsControlMycotoxin_")) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc <= -lfc_thr) %>%
  arrange(lfc)

write.xlsx(list("AFB1vsCON_up" = AFB1vsCON_up,
                "AFB1vsCON_down" = AFB1vsCON_down,
                "DONvsCON_up" = DONvsCON_up,
                "DONvsCON_down" = DONvsCON_down,
                "FB1vsCON_up" = FB1vsCON_up,
                "FB1vsCON_down" = FB1vsCON_down,
                "DONvsAFB1_up" = DONvsAFB1_up,
                "DONvsAFB1_down" = DONvsAFB1_down,
                "FB1vsAFB1_up" = FB1vsAFB1_up,
                "FB1vsAFB1_down" = FB1vsAFB1_down,
                "FB1vsDON_up" = FB1vsDON_up,
                "FB1vsDON_down" = FB1vsDON_down,
                "MYCOvsCON_up" = MYCOvsCON_up,
                "MYCOvsCON_down" = MYCOvsCON_down),
           file = file.path(outdir, "abundance.xlsx"))

