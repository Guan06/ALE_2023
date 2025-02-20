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

p1_c <- ggplot(all_sig, aes(Passage, Fold_change, group = Group)) + 
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.2, ymax=0.2, 
           alpha=0.3, fill="gold") + 
  geom_point(aes(shape = FDR_sig), color = "gray74", size = 2.4, alpha = 0.6) +
  geom_point(data = all_sig[all_sig$Compound %in% sig_lst, ],
             aes(Passage, Fold_change, shape = FDR_sig, color = Compound),
             alpha = 0.9, size = 3.4) +
  scale_shape_manual(values = c(1, 16)) +
  geom_line(aes(linetype = Concentration), color = "gray", 
            linewidth = 1.1, alpha = 0.6) + 
  geom_line(data = all_sig[all_sig$Compound %in% sig_lst, ],
            aes(Passage, Fold_change, 
                group = Group, color = Compound, linetype = Concentration),
            alpha = 0.9, linewidth = 1.1) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid")) +
  geom_hline(yintercept = 0, linetype = "twodash", 
             color = "gold", linewidth = 0.8) +
  guides(colour = guide_legend(order = 1, ncol = 1), 
         shape = guide_legend("Significance (FDR)", order = 3, ncol = 1),
         linetype = guide_legend(order = 2, ncol = 1)) +
  labs(y = "Fold change") +
  scale_color_manual(values = c("#39b87f", "#ba7e45", "#ea8783")) +
  main_theme + theme(legend.position = "right")

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

p1_d <- ggplot(lma_dat, aes(interaction(Concentration,Compound, Plate, Passage),
                            auc_l, fill = Compound)) +
  geom_boxplot(aes(color = Compound, alpha = Concentration),
               outlier.shape = NA,
               position = position_dodge2(preserve = "single")) + #,
  #   show_guide = FALSE) +
  #scale_x_discrete(labels = rep(c("0","4", "7", "10", "14", "17", "20"),
  #                              each = 8)) +
  geom_point(aes(shape = Plate, color = Compound),
             position = position_jitterdodge()) +
  main_theme +
  scale_alpha_manual(values = c("0" = 0.5,
                                "25" = 0.1,
                                "50" = 0.6)) +
  scale_color_manual(values = c("DMSO_control" = "#E69F00",
                                "Ezetimibe"= "#39b87f",
                                "Rosuvastatin" = "#b07aa1", 
                                "Simvastatin"= "#ea8783")) +
  scale_fill_manual(values = c("DMSO_control" = "#E69F00",
                               "Ezetimibe"= "#39b87f",
                               "Rosuvastatin" = "#b07aa1", 
                               "Simvastatin"= "#ea8783"), guide = "none") +
  scale_shape_manual(values = c(2, 3)) +
  labs(y = "AUC of logistic curve", 
       x = "Passage 0 -> 4 -> 7 -> 10 -> 14 -> 17 -> 20") +
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 3),
         alpha = guide_legend(order = 2)) +
  theme(legend.position = "right", axis.text.x = element_blank())

################################################################################
lma_sig <- all_sig[all_sig$Compound %in% lma_lst, ]
write.table(lma_sig, "../04.results/Figure_1d_sig.txt",
            quote = F, row.names = F, sep = "\t")

################################################################################
p1_cd <- plot_grid(p1_c, p1_d, nrow = 2, align = "v", axis = "lr",
                   rel_heights = c(1, 1),
                   labels = c("c", "d"))

ggsave("../05.figures/Figure_1_cd.pdf", p1_cd, width = 10, height = 6)
