### Make sure to run this before running Figure_S5.R
source("settings.R")
# read in all variants detected in ALE 1.0 and ALE 2.0
all1 <- readRDS("../00.data/vcf_20230829_all_variants_AF_Ratio_ALE1.rds")
all2 <- readRDS("../00.data/vcf_20230814_all_variants_ALE2.rds")
des <- read.table("../00.data/vcf_20221003_des.txt", header = T, sep = "\t")
des <- unique(des[, -2])

### Assign Plate as 2 for Metformin, as it was sequenced in several plates but 
### during experiment it was in plate 2
### Assign Plate 2 to DMSO_control that was sequenced in other plates
des[des$Compound == "Metformin", ]$Plate <- 2
des[(des$Compound == "DMSO_control" & des$Column == 9), ]$Plate <- 2
des[(des$Compound == "DMSO_control" & des$Well == "A7"), ]$Plate <- 2
write.table(des, "../04.results/20231024_vcf_des_reformat.txt", quote = F,
            sep = "\t", row.names = F)

common <- intersect(colnames(all1), colnames(all2))
all1_sub <- all1[, common]
all1_sub$Exp <- "ALE1"
all2_sub <- all2[, common]
all2_sub$Exp <- "ALE2"

################################################################################
###################### parental variants identified in pre_Figure_2.R
parental_lst <- c("1166763", "1681625", "1702756", "3356587")
###################### common variants identified in pre_Figure_2.R
com_tab <- read.table("../06.tables/Appendix_tab_01_common_variants_ALE2.txt",
                      header = T, sep =  "\t")
com_lst <- unique(com_tab$POS)
all2_sub <- all2_sub[!(all2_sub$POS %in% com_lst), ]

both <- rbind(all1_sub, all2_sub)
both <- both[!(both$POS %in% parental_lst), ]
#### 323 ALE samples (including 2 parental) + 12 ALE2 = 335 samples

rm(all1_sub, all2_sub)
gc()
###############################################################################
################################################## Number of variants AF > 0.05
### Split the variation according to the number of alternative allele
all_single <- both[both$Alt_Num == 1, ]
all_other <- both[both$Alt_Num > 1, ]

all_single_005 <- all_single[all_single$Ratio > 0.05, ]

all_other_stat <- unique(all_other %>% group_by(Sample_ID, POS, Ref) %>%
                           reframe(Ratio = sum(Alt_depth) / 
                                     (sum(Alt_depth) + Ref_depth)))

all_other_stat_005 <- all_other_stat[all_other_stat$Ratio > 0.05, ]

all_other_005 <- merge(all_other,
                       all_other_stat_005[, c("Sample_ID", "POS", "Ref")])


all_005 <- rbind(all_single_005, all_other_005)
saveRDS(all_005, "../04.results/20231023_all_005.rds")

all_stat2 <- all_005 %>% group_by(Exp, Sample_ID) %>% 
  summarise(Number_of_variants = length(POS))

all_stat2 <- as.data.frame(all_stat2)

tab <- merge(all_stat2, des)
tab$Concentration <- as.character(tab$Concentration)
tab$Plate <- as.character(tab$Plate)
tab <- tab[tab$Number_of_variant > 0, ]
write.table(tab, "../04.results/Figure_S5_variant_number_AF_005.txt",
            quote = F, sep = "\t", row.names = F)

###############################################################################
################################################## Number of variants AF > 0.5
### Split the variation according to the number of alternative allele
all_single_05 <- all_single[all_single$Ratio > 0.5, ]

all_other_stat <- unique(all_other %>% group_by(Sample_ID, POS, Ref) %>%
                           reframe(Ratio = sum(Alt_depth) / 
                                     (sum(Alt_depth) + Ref_depth)))
all_other_stat_05 <- all_other_stat[all_other_stat$Ratio > 0.5, ]

all_other_05 <- merge(all_other,
                      all_other_stat_05[, c("Sample_ID", "POS", "Ref")])


all_05 <- rbind(all_single_05, all_other_05)
saveRDS(all_05, "../04.results/20231023_all_05.rds")

all_stat2 <- all_05 %>% group_by(Exp, Sample_ID) %>% 
  summarise(Number_of_variants = length(POS))

all_stat2 <- as.data.frame(all_stat2)

tab <- merge(all_stat2, des)
tab$Concentration <- as.character(tab$Concentration)
tab$Plate <- as.character(tab$Plate)
#tab <- tab[tab$Number_of_variant > 0, ]
write.table(tab, "../04.results/Figure_S5_variant_number_AF_05.txt",
            quote = F, sep = "\t", row.names = F)