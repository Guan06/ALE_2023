################ make sure to run pre_Figure_1.R before and have the output 
################ generated as ../04.results/Figure_1a_growth_curve.txt
## source("pre_Figure_1.R")
source("plot_settings.R")
source("function_get_sig.R")

##############################################################################
################ Compare each compound vs control
all_dat <- read.table("../04.results/Figure_1_cd_growth_curve.txt",
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
ale1 <- read.table("../04.results/Figure_S3_growth_curve_24h_OD.txt",
                   header = T, sep = "\t")
####################Read in OD and plate layout for ALE 2.0
ale2_raw <- read.table("../00.data/growth_curve_24h/ALE2_20230125_P0_P5_P10_P15_P20.txt",
                   header = T, sep = "\t")
well2 <- read.table("../00.data/growth_plate_layout_20230125.txt",
                    header = T, sep = "\t")
well2 <- well2[, c("Well", "Compound", "Concentration", "Solvent")]


ale2 <- ale2_raw %>% pivot_longer(!c(Date, Passage, Time),
                                  names_to = "Well",
                                  values_to = "OD") 
ale2 <- merge(ale2, well2)
ale2$ID <- paste0("ALE2_", ale2$Well)

###### Merge and plot PFOA, PFNA
ale12 <- rbind(ale1[, c("ID", "Time", "OD", "Compound",
                        "Passage", "Concentration")],
               ale2[, c("ID", "Time", "OD", "Compound",
                        "Passage", "Concentration")])
ale12 <- ale12[ale12$Time < 24, ]
ale12 <- ale12[ale12$Compound %in% c("PFNA", "PFOA", "DMSO_control"),]

## remove the data with OD > 1 -> most possibly a mistake by machine
ale12 <- ale12[ale12$OD < 1, ]
ale12$Concentration <- as.character(ale12$Concentration)

### Plot the passages that overlapped between ALE1 and ALE2 in main figures
### Passage 0, 10, 20
ale12_main <- ale12[ale12$Passage %in% c(0, 10, 20), ]

p1_d <- ggplot(ale12_main[ale12_main$Compound == "DMSO_control", ],
               aes(Time, OD, group = ID)) +
  geom_line(aes(color = Compound), alpha = 0.6) +
  geom_point(aes(color = Compound, shape = Concentration),
            alpha = 0.6, size = 2) +
  geom_line(data = ale12_main[ale12_main$Compound != "DMSO_control", ],
            aes(color = Compound), alpha = 0.8) +
  geom_point(data = ale12_main[ale12_main$Compound != "DMSO_control", ],
             aes(color = Compound, shape = Concentration),
             alpha = 0.8, size = 2) +
  scale_shape_manual(values = c("0" = 1, "25" = 10, "50" = 19,
                                "250" = 2, "500" = 17)) +
  facet_wrap(~Passage, scales = "fixed", nrow = 1) +
  main_theme +
  labs(x = "Time / h") +
  scale_color_manual(values = c("DMSO_control" = "#E69F00",
                                 "PFNA" = "lightblue4",
                                 "PFOA" = "yellow4"))  +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, 
          size = rel(1.1), face = "bold")) 

p1_cd <- plot_grid(p1_c, p1_d, nrow = 2, align = "v", axis = "lr",
                   rel_heights = c(1, 1.1),
                   labels = c("c", "d"))

ggsave("../05.figures/Figure_1_cd.pdf", p1_cd, width = 10, height = 6)