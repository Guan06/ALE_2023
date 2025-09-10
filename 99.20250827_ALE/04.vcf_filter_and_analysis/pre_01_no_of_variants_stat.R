setwd("~/Desktop/ALE_2023/99.20250827_ALE/04.vcf_filter_and_analysis/")
library(tidyverse)

v2023 <- readRDS("../../00.data/vcf_20230829_all_variants_AF_Ratio_ALE1.rds")
#v2025 <- readRDS("../03.merge_vcf/20250907_all_variants_redo.rds")
v2025 <- readRDS("../03.merge_vcf/20250910_all_variants_redo_qual.rds")

v2023_stat <- v2023 %>% group_by(Sample_ID) %>% summarise(No_of_var_2023 = n())
v2025_stat <- v2025 %>% group_by(Sample_ID) %>% summarise(No_of_var_2025 = n())

stat <- merge(v2023_stat, v2025_stat)
stat[stat$No_of_var_2023 != stat$No_of_var_2025, ]$Sample_ID
# Sample_ID No_of_var_2023 No_of_var_2025
# Plate1E4            360            645

des <- read.delim("../03.merge_vcf/20220916_design.txt", header = T, sep = "\t")
des[des$Sample_ID == "Plate1E4", ]

v2023_1E4 <- v2023[v2023$Sample_ID == "Plate1E4", ]
rm(v2023)
