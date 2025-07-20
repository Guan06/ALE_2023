################ make sure to run pre_Figure_1.R before and have the output 
################ generated as ../04.results/Figure_1cde_growth_curve.txt
## source("pre_Figure_1.R")
source("plot_settings.R")
source("function_get_sig.R")

##############################################################################
################ Compare each compound vs control
all_dat <- read.table("../04.results/Figure_1_cde_growth_curve.txt",
                      header = T, sep = "\t")
all_sig <- c()
all_comp_lst <- unique(all_dat$Compound)
for (c in all_comp_lst) {
  if (c %in% c("DMSO_control", "Water_control", "Negative_control", 
               "Cumic_alcohol")) {next}

  this_c <- all_dat[all_dat$Compound == c, ]
  this_plate <- unique(this_c$Plate)
  if (length(unique(this_c$Solvent)) > 1) {
    this_solvent <- "DMSO_control"
  } else if (grepl(unique(this_c$Solvent), "Water")) {
    this_solvent <- "Water_control"
  } else {
    this_solvent <- "DMSO_control"
  }
  #this_solvent <- paste0(unique(this_c$Solvent), "_control")
  this_s <- all_dat[all_dat$Compound == this_solvent, ]
  
  this <- rbind(this_c, this_s)
  
  ################# significant test and FC
  ## Only keep the DMSO control from the same plate
  this <- this[this$Plate %in% this_plate, ]
  
  this$Group <- paste0(this$Passage, "_", this$Compound,
                       "_", this$Concentration)
  
  sig <- get_sig(this, "Group", "auc_l")
  
  sig$Compound <- c
  all_sig <- rbind(all_sig, sig)
}

all_sig$Fold_change <- as.numeric(all_sig$Fold_change)
all_sig_hit <- all_sig[(abs(all_sig$Fold_change) >= 0.20) &
                      (all_sig$FDR < 0.05), ]
write.table(all_sig, "../04.results/Figure_1c_all_gc_cmp.txt", quote = F,
            sep = "\t", row.names = F)
write.table(all_sig_hit, "../04.results/Figure_1c_hit_gc_cmp.txt", quote = F,
            sep = "\t", row.names = F)
###############################################################################

all_sig$Group <- paste0(all_sig$Compound, "_", all_sig$Concentration)
all_sig$Passage <- as.integer(all_sig$Passage1)
all_sig$Concentration <- as.character(all_sig$Concentration)
################################################################################
### not presenting ultra-high concentration PFAS here
all_sig <- all_sig[all_sig$Concentration %in% c("25", "50"), ]
all_sig_hit <- all_sig[(abs(all_sig$Fold_change) >= 0.20) &
                         (all_sig$FDR < 0.05), ]
sig_lst <- unique(all_sig_hit$Compound)

library(colorblindcheck)
#palette_check(c("#20b2aa", "#7a67ee", "#FA8072"), plot = TRUE)
palette_check(c("#ffc0cb", "#98dd94", "#6495ed"), plot = TRUE)

p1_c <- ggplot(all_sig, aes(Passage, Fold_change, group = Group)) + 
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.2, ymax=0.2, 
           alpha=0.15, fill="gold") + 
  geom_point(aes(shape = FDR_sig), color = "gray74", size = 2.4, alpha = 0.6) +
  scale_shape_manual(values = c(1, 16)) +
  geom_line(aes(linetype = Concentration), color = "gray", 
            linewidth = 1.1, alpha = 0.6) + 
  geom_line(data = all_sig[all_sig$Compound %in% sig_lst, ],
            aes(Passage, Fold_change, 
                group = Group, color = Compound, linetype = Concentration),
            alpha = 0.9, linewidth = 1.1) +
  geom_point(data = all_sig[all_sig$Compound %in% sig_lst, ],
             aes(Passage, Fold_change, shape = FDR_sig, color = Compound),
             alpha = 0.9, size = 3.4) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid")) +
  geom_hline(yintercept = 0, linetype = "twodash", 
             color = "gold", linewidth = 0.8) +
  guides(colour = guide_legend(order = 1, ncol = 1), 
         shape = guide_legend("Significance (FDR)", order = 3, ncol = 1),
         linetype = guide_legend(order = 2, ncol = 1)) +
  labs(y = "Fold change") +
  scale_color_manual(values = c("Ezetimibe" = "#ffc0cb",
                                "Loperamide" = "#98dd94",
                                "Simvastatin" = "#6495ed"))+
  #  "#FA8072", "#20b2aa", "#7a67ee")) + 
  #"#E69F00", "#8BC34A", "#FFC107",
  #"#FFA07A", "#34A85A", "#FFD700", "#964B00"
  ## old color for loperamide is #ba7e45
  main_theme + theme(legend.position = "right")

p1_c
##############################################################################

##############################################################################
############### Compare lipid modification agents
all_dat <- read.table("../04.results/Figure_1_cde_growth_curve.txt",
                      header = T, sep = "\t")
lma_lst <- c("Ezetimibe", "Rosuvastatin", "Simvastatin", "DMSO_control")
lma_dat <- all_dat[all_dat$Compound %in% lma_lst, ]
## only take DMSO_control from Plate_2 and Plate_3
lma_dat <- lma_dat[lma_dat$Plate %in% c("Plate_2", "Plate_3"), ]

lma_dat$Passage <- as.integer(lma_dat$Passage)
lma_dat$Concentration <- as.character(lma_dat$Concentration)

lma_dat_DMSO <- lma_dat[lma_dat$Compound == "DMSO_control", ]
lma_dat_DMSO_mean <- lma_dat_DMSO %>% group_by(Passage, Plate) %>% 
  summarise(auc_l_DMSO = mean(auc_l))

lma_dat_cmp <- lma_dat[lma_dat$Compound != "DMSO_control", ]

dat_1d <- merge(lma_dat_cmp, lma_dat_DMSO_mean, by = c('Plate', "Passage"))

palette_check(c("#ccccff", "#6495ed", "#ffc0cb", "#E69F00"), plot = TRUE)


p1_d <- ggplot(lma_dat, aes(interaction(Concentration,Compound, Plate, Passage),
                            auc_l, fill = Compound)) +
  geom_boxplot(aes(color = Compound, alpha = Concentration),
               outlier.shape = NA,
               position = position_dodge2(preserve = "single")) + #,
  #   show_guide = FALSE) +
  #scale_x_discrete(labels = rep(c("0","4", "7", "10", "14", "17", "20"),
  #                              each = 8)) +
  geom_point(aes(shape = Plate, color = Compound),
             position = position_jitterdodge(), alpha = 1) +
  main_theme +
  scale_alpha_manual(values = c("0" = 0,
                                "25" = 0,
                                "50" = 0.8)) +
  # scale_color_manual(values = c("DMSO_control" = "#E69F00",
  #                               "Ezetimibe"= "#20b2aa",
  #                               "Rosuvastatin" = "#9f1717", 
  #                               "Simvastatin"= "#7a67ee")) +
  # 
  scale_color_manual(values = c("DMSO_control" =  "#E69F00",
                                "Ezetimibe" = "#ffc0cb",
                                "Rosuvastatin" = "#ccccff",
                                "Simvastatin" = "#6495ed")) +
  scale_fill_manual(values = c("DMSO_control" =  "#E69F00",
                                "Ezetimibe" = "#ffc0cb",
                                "Rosuvastatin" = "#ccccff",
                                "Simvastatin" = "#6495ed"), guide = "none") +
  
  # scale_fill_manual(values = c("DMSO_control" = "#E69F00",
  #                              "Ezetimibe"= "#20b2aa",
  #                              "Rosuvastatin" = "#9f1717", 
  #                              "Simvastatin"= "#7a67ee"), guide = "none") +
  scale_shape_manual(values = c(2, 3)) +
  labs(y = "AUC of logistic curve", 
       x = "Passage 0 -> 4 -> 7 -> 10 -> 14 -> 17 -> 20") +
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 3),
         alpha = guide_legend(order = 2)) +
  theme(legend.position = "right", axis.text.x = element_blank())

p1_d
################################################################################
lma_sig <- all_sig[all_sig$Compound %in% lma_lst, ]
write.table(lma_sig, "../04.results/Figure_1d_sig.txt",
            quote = F, row.names = F, sep = "\t")

################################################################################
p1_cd <- plot_grid(p1_c, p1_d, nrow = 2, align = "v", axis = "lr",
                   rel_heights = c(1, 1),
                   labels = c("c", "d"))
p1_cd

################################################################################
## Plate 1 and Plate 2 contains mGAM + DMSO
## Plate 3 and Plate 4 contains mGAM + DMSO + 50 µM Loperamide

this_dat <- read.table("../00.data/growth_curve_24h/tidy_20250514_Buni_Loperamide.txt",
                       header = T, sep = "\t")
well <- read.table("../00.data/growth_plate_layout_20250515.txt", 
                   header = T, sep = "\t")

well_order <- c("Time", well$Well[1:96])

plate_lst <- unique(this_dat$Plate)
plate_df <- data.frame(Plate = c("Plate_1", "Plate_2", "Plate_3", "Plate_4"),
                       Condition = c("DMSO", "DMSO",
                                     "Loperamide", "Loperamide"))
library(growthcurver)
all_gc <- c()
for (pl in plate_lst){
  this_pl <- this_dat[this_dat$Plate == pl, ]
  this_pl <- this_pl[, -1]
  this_pl <- this_pl[, match(well_order, colnames(this_pl))]
  
  if (nrow(this_pl) >24) {this_pl <- this_pl[1:24, ]}
  
  this_gc <- SummarizeGrowthByPlate(this_pl, plot_fit = F)
  this_gc$Plate <- pl
  all_gc <- rbind(all_gc, this_gc)
}

colnames(all_gc)[1] <- "Well"
all_dat <- merge(all_gc, well)
all_dat_noNC <- all_dat[all_dat$Isolate != "NC", ] 
all_dat_noNC <- merge(all_dat_noNC, plate_df)

### Filter out Plate 1 G6 well as we saw some unsolvable material in that well
### Could be either paper or plastic pieces
all_dat_noNC <- all_dat_noNC[!(all_dat_noNC$Plate == "Plate_1" &
                                 all_dat_noNC$Well == "G6"), ]
all_dat_noNC$Row_num <- as.factor(all_dat_noNC$Row_num)
all_dat_noNC$Group <- factor(all_dat_noNC$Group,
                             levels = c("Parental", "DMSO_evo", "Lop50_evo"),
                             ordered = TRUE)
p1e_right <- ggplot(all_dat_noNC, 
                    aes(interaction(Row_num, Group), auc_l)) +
  # geom_boxplot(aes(color = Group),
  #              outlier.shape = NA,
  #              position = position_dodge2(preserve = "single")) + 
  geom_boxplot(aes(color = Group), outlier.shape = NA) + 
  geom_point(aes(color = Group, shape = Plate),
             position = position_jitterdodge(), alpha = 0.8) +
  facet_wrap(~Condition) +
  scale_color_manual(values = c("Lop50_evo" = "#98dd94",
                                "Parental" = "#808060",
                                "DMSO_evo" = "#E69F00")) +
  #scale_shape_manual(values = 1:8) +
  main_theme +
  theme(axis.text.x = element_text(colour = "black", angle = 90,
                                   size = 6, hjust = 1),
        legend.position = "NA") +
  theme(axis.text.x = element_blank()) +
  scale_shape_manual(values = c(4, 1, 4, 1)) +
  labs(y = "AUC of logistic curve", 
       x = "Isolates 1 to 8 of parental strain, DMSO- and Loperamide-evolved strain") 

p1e_right

################################################################################
### Significant test for panel e
sig <- c()

for (c in unique(all_dat_noNC$Condition)) {
  this <- all_dat_noNC[all_dat_noNC$Condition == c,]
  this_sig <- get_sig_iso(this, "Group", "auc_l")
  this_sig$Condition <- c
  sig <- rbind(sig, this_sig)
}

sig$FDR <- p.adjust(sig$Significance, method = "fdr")
sig$FDR_sig <- ifelse(sig$FDR < 0.05, "Sig", "Non-sig")

write.table(sig, "../04.results/Figure_1e_sig.txt",
            quote = F, row.names = F, sep = "\t")
################################################################################

library(agricolae)
for (c in unique(all_dat_noNC$Condition)) {
  this <- all_dat_noNC[all_dat_noNC$Condition == c,]
  print(c)
  print(kruskal(this$auc_l, this$Group, p.adj = "fdr"))
}

##### Plot the growth curve
od <- this_dat %>% pivot_longer(!c("Plate", "Time"), 
                                names_to = "Well", values_to = "OD")
od <- merge(od, well)
od <- od[od$Group != "NC", ]

### Filter out Plate 1 G6 well as we saw some unsolvable material in that well
### Could be either paper or plastic pieces
od <- od[!(od$Plate == "Plate_1" & od$Well == "G6"), ]
od <- merge(od, plate_df)
od$ID2 <- paste0(od$Plate, "_", od$ID) 

p1e_left <- ggplot(od, aes(Time, OD, group = ID2)) +
  geom_line(aes(color = Group), alpha = 0.7) +
  geom_point(aes(color = Group, shape = Plate),
             alpha = 0.5, size = 1.2) +
  facet_wrap(~Condition, scales = "fixed", nrow = 1) +
  main_theme +
  scale_color_manual(values = c("Lop50_evo" = "#98dd94",
                                "Parental" = "#808060",
                                "DMSO_evo" = "#E69F00")) +
  scale_shape_manual(values = c(4, 1, 4, 1))  

p1e_left

p1e <- plot_grid(p1e_left, p1e_right, nrow = 1, labels = c("e", ""),
                 align = "h", rel_widths = c(2, 1))
p1e

### Add flow cytometry result to panel f
fc <- read.table("../00.data/20250611_FACS.txt", header = T, sep = "\t")
fc$ID3 <- paste0(fc$Plate, fc$Well)

fc_od <- od[(od$Plate %in% fc$Plate) & od$Time == 24, ]
fc_od$ID3 <- paste0(fc_od$Plate, fc_od$Well)

fc_od <- merge(fc, fc_od)
p1f <- ggplot(fc_od, aes(OD, cell.mL_dilution_factor)) +
  geom_smooth(method="lm", col="gray82", se = FALSE, alpha = 0.8) + 
  geom_point(aes(color = Group, shape = Condition), size = 2) +
  scale_color_manual(values = c("Lop50_evo" = "#98dd94",
                                "Parental" = "#808060",
                                "DMSO_evo" = "#E69F00")) +
  main_theme +
  scale_shape_manual(values = c(1, 16)) +
  # theme(axis.text.x = element_text(colour = "black", angle = 90,
  #                                  size = 6, hjust = 1),
  #       legend.position = "NA") +
  #theme(axis.text.x = element_blank()) +
  labs(y = "Cells / mL ", 
       x = "OD") 
p1f
cor.test(fc_od$OD, fc_od$cell.mL_dilution_factor, method = "spearman")

env_lst <- unique(fc$Condition)
for (c in env_lst) {
  this_fc <- fc[fc$Condition == c, ]
  
  
}
################################################################################
### Significant test for panel f
# sig <- c()
# 
# for (c in unique(fc$Condition)) {
#   this <- fc[fc$Condition == c,]
#   this_sig <- get_sig_iso(this, "Group", "cell.mL_dilution_factor")
#   this_sig$Condition <- c
#   sig <- rbind(sig, this_sig)
# }
# 
# sig$FDR <- p.adjust(sig$Significance, method = "fdr")
# sig$FDR_sig <- ifelse(sig$FDR < 0.05, "Sig", "Non-sig")
# 
# write.table(sig, "../04.results/Figure_1f_sig.txt",
#             quote = F, row.names = F, sep = "\t")
################################################################################

plot_grid(p1e, p1f, rel_widths = c(2, 1))

p1_cde <- plot_grid(p1_cd, p1e, nrow = 2, align = "v", axis = "l",
                    rel_heights = c(1.9, 1))

ggsave("../05.figures/Figure_1_cde.pdf", p1_cde, width = 10, height = 8)
