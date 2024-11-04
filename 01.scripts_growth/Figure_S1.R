source("plot_settings.R")

library(rcdk)
library(dplyr)
library(tidyr)
library(pheatmap)

lst <- read.delim("../00.data/20240821_chem_SMILEs.txt", header = T, sep = "\t")
## generate morgen finger prints from SMILE

nlst <- nrow(lst)
all <- c()
for (i in 1:nlst) {
  for (j in 1:nlst){
    moles <- parse.smiles(lst[c(i, j), 2])
    fps <- lapply(moles, get.fingerprint, type = "extended",
                  fp.mode = "bit", size = 1024)
    
    
    fp1 <- fps[[1]]
    fp2 <- fps[[2]]
    
    inter_count <- length(intersect(fp1@bits, fp2@bits))
    union_count <- length(union(fp1@bits, fp2@bits))
    
    tanimoto <- inter_count / union_count
    all <- rbind(all, c(lst[i, 1], lst[j, 1], tanimoto))
  }
}

colnames(all) <- c("Compound_1", "Compound_2", "Tanimoto_coefficient")

all <- as.data.frame(all)
all$Tanimoto_coefficient <- as.numeric(all$Tanimoto_coefficient)

##########################
### Plot heatmap of compounds with category information

all_mat <- all %>% 
  pivot_wider(names_from = Compound_2, values_from = Tanimoto_coefficient)

tmp <- all_mat$Compound_1

all_mat <- as.matrix(all_mat[, -1])
rownames(all_mat) <- tmp


meta <- read.table("../00.data/compound_category.txt",
                   header = T, sep = "\t")
meta <- rbind(meta, c("DMSO", "Solvent"))

meta2 <- data.frame(Category = factor(x = meta$Category, 
                                      levels = unique(meta$Category)))
rownames(meta2) <- meta$Compound

## Define colors
mat_colors <- list(Category = palette_OkabeIto[1:length(unique(meta2$Category))])
names(mat_colors$Category) <- unique(meta2$Category)

# all_mat[all_mat < 0.85] <- 0
# Only PFNA/PFOA and Cypermethrin/Permethrin pair have > 0.85 Tanimoto coefficient 

#all_mat[lower.tri(all_mat)] <- NA

p <- pheatmap(all_mat, 
              annotation_col = meta2, annotation_colors = mat_colors,
              annotation_row = meta2,
              treeheight_row = 0, treeheight_col = 0,
              na_col="white",
              border_color ="gray")

ggsave("../05.figures/Figure_S1a.pdf", p, height = 8, width = 10.7)
