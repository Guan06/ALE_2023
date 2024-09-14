source("settings.R")

# read in all variants detected in ALE 1.0 and ALE 2.0 with AF > 0.05
all <- readRDS("../03.results/20231023_all_005.rds")
all[all$Compound == "DMSO", ]$Compound <- "DMSO_control"

###############################################################################
########## Get thevariant with TYPE snp or snp,snp
all_snp <- all[all$TYPE %in% c("snp", "snp", "snp"), ]
# 74.42%
all_pos <- unique(all_snp[, c("POS", "REF")])
write.table(all_pos, "../03.results/Figure_S8_mut_pos.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)
### --> export to process with perl script get_tri_nt.pl
############# perl get_tri_nt.pl ../03.results/Figure_S8_mut_pos.txt >
###  ../03.results/Figure_S8_mut_tri_nt.txt

###############################################################################
tri <- read.table("../03.results/Figure_S8_mut_tri_nt.txt", 
                  header = T, sep = "\t")
colnames(tri)[1] <- "POS"

### Get the mutational profile for each sample
all_snp <- merge(all_snp, tri)
all_snp$Ref_nt <- ifelse(all_snp$REF %in% c("C", "T"), all_snp$REF, 
                          chartr("AG", "TC", all_snp$REF))


all_snp$Alt_nt <- ifelse(all_snp$REF %in% c("C", "T"), all_snp$ALT,
                          chartr("ATCG", "TAGC", all_snp$ALT))

all_snp$Tri_nt <- ifelse(all_snp$REF %in% c("C", "T"), all_snp$Trinucleotide,
                          all_snp$RC_trinucleotide)

all_snp$ID <- paste0(all_snp$Ref_nt, "_", all_snp$Alt_nt, ".",
                      all_snp$Tri_nt)

id_to_mut <- read.table("../00.data/20221127_ID_to_mutation_type.txt", 
                        header = T, sep = "\t")
all_snp <- merge(all_snp, id_to_mut)

write.table(all_snp, "../03.results/Figure_S8_all_snp.txt", 
            quote = F, sep = "\t", row.names = F)
###############################################################################
########## Get Alt_Num == 1 for SNPs

#all %>% group_by(Alt_Num) %>% summarise(n())
## N = 1, 9317 (90.2%)
## N = 2, 900 (8.66)
## N = 3, 117 (1.13%)

#all_1 <- all[all$Alt_Num == 1, ]
# complex   131
# del       497
# ins      1008
# mnp         4
# snp      7731 (82.5%)
#all_1_snp <- all_1[all_1$TYPE == "snp", ]
#pos1 <- unique(all_1_snp[, c("POS", "REF")])

###############################################################################
