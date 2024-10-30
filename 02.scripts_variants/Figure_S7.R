## copy number of SusC across genomes
## results in RDS
## /rds-d6/user/rg684/hpc-work/ALE_20230811_batch2/03.MAGs/00.TonB
library(forcats)
source("settings.R")

meta <- read.table("../00.data/Bacteroides_uniformis_metadata.tsv", header = T,
                   sep = "\t")
#############################################################
## filter out 10 genomes without country/continent information
meta <- meta[!is.na(meta$Continent), ]
meta <- meta[meta$Continent != "not provided", ]
# 5990 genomes left

meta_count <- as.data.frame(meta %>% group_by(Genome_type, Country) %>% 
                              summarise(Number = length(Genome)))

## filter out countries with less than 5 genomes
lst <- meta_count[meta_count$Number < 5, ]$Country
meta <- meta[!(meta$Country %in% lst), ]
# 5984, filtered out "Bangladesh" "Madagascar" "Singapore" 

## Add copy number of SusC
sus <- read.table("../00.data/copy_number_SusC.txt", header = T, sep = "\t")
meta <- merge(meta, sus)

meta$Country <- fct_reorder(meta$Country, .x = meta$Length,
                            .fun = median, .desc = T)

meta$Quality <- ifelse((meta$Completeness > 90 & meta$Contamination < 5), 
                       "High", 
                       ifelse((meta$Completeness >50 & meta$Contamination < 10),
                              "Medium", low))
## One MAG has size 7700188, filter out as completeness only 58.62
meta <- meta[meta$Length < 7700188, ]

p1 <- ggplot(meta, aes(Country, Length)) +
#  geom_violin(aes(color = Continent)) +
  geom_boxplot(width=0.2, color = "gray38", fill = "NA", outlier.shape = NA) +
  geom_violin(aes(color = Continent, fill = Continent), alpha = 0.5) +
  geom_jitter(aes(shape = Quality), color = "gray58", 
              width = 0.2, size = 1, alpha = 0.6) +
  scale_shape_manual(values = c(1, 3)) +
  scale_color_manual(values = palette_OkabeIto[c(1:3, 7)]) +
  scale_fill_manual(values = palette_OkabeIto[c(1:3, 7)]) +
  geom_hline(yintercept = 4683603, linetype = "dashed", color = "gray") +
  #facet_wrap(~Genome_type, ncol = 1, scales = "free_y") +
  facet_grid(~Genome_type, scales = "free", space = "free") +
  main_theme + ylab("Genome size") + xlab("") + 
  theme(axis.text.x = element_blank())

p2 <- ggplot(meta, aes(Country, Copy_number)) +
  geom_violin(aes(color = Continent, fill = Continent), alpha = 0.5) +
  geom_jitter(aes(shape = Quality), color = "gray58", 
              width = 0.2, size = 1, alpha = 0.6) +
  #geom_boxplot(width=0.2, color = "gray38", outlier.shape = NA) +
  geom_hline(yintercept = 61, linetype = "dashed", color = "gray") +
  scale_color_manual(values = palette_OkabeIto[c(1:3, 7)]) +
  scale_fill_manual(values = palette_OkabeIto[c(1:3, 7)]) +
  scale_shape_manual(values = c(1, 3)) +
  #facet_wrap(~Genome_type, ncol = 1, scales = "free_y") +
  facet_grid(~Genome_type, scales = "free", space = "free") +
  main_theme + ylab("Copy number of susC") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#p2
p <- plot_grid(p1, p2, nrow = 2, align = "hv", axis = "l", 
               rel_heights = c(1, 1))

ggsave("../04.figures/Figure_S7.pdf", width = 10, height = 8)

summary(aov(Copy_number ~ Length * Contamination * Genome_type * Country,
            meta))

ggplot(meta, aes(Length, Copy_number, color = Continent)) +
  geom_point(shape = 1) + 
  scale_color_manual(values = palette_OkabeIto[c(1:3, 7)]) +
  geom_smooth(method = "lm")
