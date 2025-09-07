#!/usr/local/Cluster-Apps/R/R.4.0.3/bin/Rscript

samples <- list.dirs("../02.call_variant_snippy/mincov10_C3_F001/", full.names = F, recursive = F)

all <- c()

## These are the samples that have too low sequencing depth; including
## 1 water-control, 2 DMSO-control, and 4 xenobiotic-evolved pop

skip_lst <- c("Plate3B12", "Plate3H7", "Plate3G12",
              "Plate3H11", "Plate3H12", "Plate4D6", "Plate3H6")

for (s in samples){
    if (s %in% skip_lst) {next}
    this_file <- paste0("../02.call_variant_snippy/mincov10_C3_F001/", s, "/snps.tab")
    print(this_file)

    this_tab <- read.table(this_file, header = T, sep = "\t", quote = "", fill = FALSE)

    this_tab$Sample_ID <- s

    all <- rbind(this_tab, all)
}

des <- read.table("20220916_design.txt", header = T, sep = "\t")
des <- unique(des[, -1])

all <- merge(all, des)

## remove columns that are not useful
all <- subset(all, select = -c(CHROM, GENE, LOCUS_TAG))
##write.table(all, "20250907_all_variants_redo.txt", quote = F, sep = "\t", row.names = F)

all[!grepl("CDS", all$FTYPE), ]$FTYPE <- "Non-coding"
all[!grepl("CDS", all$FTYPE), ]$EFFECT <- "Non-coding_region_variant"

library(stringr)
all$Ref <- unlist(str_split_fixed(all$EVIDENCE, " ", 2))[, 2]
all$Ref_depth <- as.numeric(unlist(str_split_fixed(all$Ref, ":", 2))[, 2])

all$Concentration <- as.character(all$Concentration)

### separate the alt allele by 1, 2, or 3 altenative alleles
all$Alt <- unlist(str_split_fixed(all$EVIDENCE, " ", 2))[, 1]

## If more than one alternative allele exhists, there will be "," in column "Alt"
all_af1 <- all[!grepl(",", all$EVIDENCE), ]
all_af1$Alt_depth <- as.numeric(unlist(str_split_fixed(all_af1$Alt, ":", 2))[, 2])
all_af1$Alt_Num <- 1
all_af1$Ratio <- all_af1$Alt_depth / (all_af1$Ref_depth + all_af1$Alt_depth)
# dim(all_af1) 483166 * 26

all_af2_af3 <- all[grepl(",", all$EVIDENCE), ]
all_af2_af3$Alt_Num <- lengths(str_split(all_af2_af3$ALT, ","))

all_af2 <- all_af2_af3[all_af2_af3$Alt_Num == 2, ]
all_af3 <- all_af2_af3[all_af2_af3$Alt_Num == 3, ]

################## Reformat POS with 2 Alt
all_af2_cp1 <- all_af2_cp2 <- all_af2
all_af2_cp1$ALT <- unlist(str_split_fixed(all_af2_cp1$ALT, ",", 2))[, 1]
all_af2_cp1$Alt_depth <- unlist(str_split_fixed(all_af2_cp1$Alt, ":", 2))[, 2]
all_af2_cp1$Alt_depth <- as.numeric(unlist(str_split_fixed(all_af2_cp1$Alt_depth,
                                                           ",", 2))[, 1])

all_af2_cp2$ALT <- unlist(str_split_fixed(all_af2_cp2$ALT, ",", 2))[, 2]
all_af2_cp2$Alt_depth <- unlist(str_split_fixed(all_af2_cp2$Alt, ":", 2))[, 2]
all_af2_cp2$Alt_depth <- as.numeric(unlist(str_split_fixed(all_af2_cp2$Alt_depth,
                                                           ",", 2))[, 2])

all_af2 <- rbind(all_af2_cp1, all_af2_cp2)
# nrow(all_af2) 1325 * 2 = 2650
all_af2$Ratio <- all_af2$Alt_depth / (all_af2$Ref_depth + all_af2$Alt_depth)

################## Reformat POS with 3 Alt
all_af3_cp1 <- all_af3_cp2 <- all_af3_cp3 <- all_af3
all_af3_cp1$ALT <- unlist(str_split_fixed(all_af3_cp1$ALT, ",", 3))[, 1]
all_af3_cp1$Alt_depth <- unlist(str_split_fixed(all_af3_cp1$Alt, ":", 2))[, 2]
all_af3_cp1$Alt_depth <- as.numeric(unlist(str_split_fixed(all_af3_cp1$Alt_depth,
                                                           ",", 3))[, 1])

all_af3_cp2$ALT <- unlist(str_split_fixed(all_af3_cp2$ALT, ",", 3))[, 2]
all_af3_cp2$Alt_depth <- unlist(str_split_fixed(all_af3_cp2$Alt, ":", 2))[, 2]
all_af3_cp2$Alt_depth <- as.numeric(unlist(str_split_fixed(all_af3_cp2$Alt_depth,
                                                           ",", 3))[, 2])

all_af3_cp3$ALT <- unlist(str_split_fixed(all_af3_cp3$ALT, ",", 3))[, 3]
all_af3_cp3$Alt_depth <- unlist(str_split_fixed(all_af3_cp3$Alt, ":", 2))[, 2]
all_af3_cp3$Alt_depth <- as.numeric(unlist(str_split_fixed(all_af3_cp3$Alt_depth,
                                                           ",", 3))[, 3])

all_af3 <- rbind(all_af3_cp1, all_af3_cp2, all_af3_cp3)
all_af3$Ratio <- all_af3$Alt_depth / (all_af3$Ref_depth + all_af3$Alt_depth)


##################### Merge all together
all <- rbind(all_af1, all_af2, all_af3)
all$Concentration <- as.character(all$Concentration)

## Add the effect and group them
all$Effect <- unlist(str_split_fixed(all$EFFECT, " ", 2))[, 1]
##saveRDS(all, "20230807_all_variants_addAF_updateRatio.rds")
saveRDS(all, "20250907_all_variants_redo.rds")
