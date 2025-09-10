setwd("~/Desktop/ALE_2023/99.20250827_ALE/04.vcf_filter_and_analysis/")
library(tidyverse)
library(ggplot2)

v2025 <- readRDS("../03.merge_vcf/20250907_all_variants_redo.rds")
v2025$Evolved <- ifelse(v2025$Compound == "NT5002", "Parental", "Evolved")

ggplot(v2025, aes(log(Ratio))) +
  geom_histogram(bins = 100) +
  theme_bw() + labs(x = "Allele Frequency (log-transformed)",
                    y = "Number of variants") +
  facet_wrap(~Evolved, scales = "free") +
  geom_vline(xintercept = log(0.05), color = "salmon", linetype = "dashed")+
  geom_vline(xintercept = log(0.5), color = "gold", linetype = "dashed")
