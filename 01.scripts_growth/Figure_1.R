################ make sure to run pre_Figure_1.R before and have the output 
################ generated as ../03.results/Figure_1a_growth_curve.txt
## source("pre_Figure_1.R")
source("plot_settings.R")
source("function_get_sig.R")

##############################################################################
################ Compare each compound vs control
all_dat <- read.table("../03.results/Figure_1a_growth_curve.txt",
                      header = T, sep = "\t")
all_sig <- c()
all_comp_lst <- unique(all_dat$Compound)
for (c in all_comp_lst) {
  if (c %in% c("DMSO_control", "Water_control", "Negative_control", 
               "Cumic_alcohol")) {next}

  this_c <- all_dat[all_dat$Compound == c, ]
  this_plate <- unique(this_c$Plate)
  this_solvent <- paste0(unique(this_c$Solvent), "_control")
  this_s <- all_dat[all_dat$Compound == this_solvent, ]
  
  this <- rbind(this_c, this_s)
  
  ################# significant test and FC
  ## Only keep the DMSO control from the same plate
  this <- this[this$Plate == this_plate, ]
  
  this$Group <- paste0(this$Passage, "_", this$Compound,
                       "_", this$Concentration)
  
  sig <- get_sig(this, "Group", "auc_l")
  
  sig$Compound <- c
  all_sig <- rbind(all_sig, sig)
}

all_sig$Fold_change <- as.numeric(all_sig$Fold_change)
all_sig_hit <- all_sig[(abs(all_sig$Fold_change) >= 0.20) &
                      (all_sig$FDR < 0.05), ]
write.table(all_sig, "../03.results/Figure_1a_all_gc_cmp.txt", quote = F,
            sep = "\t", row.names = F)
write.table(all_sig_hit, "../03.results/Figure_1a_hit_gc_cmp.txt", quote = F,
            sep = "\t", row.names = F)
###############################################################################
###############################################################################
### 
sig_lst <- unique(all_sig_hit$Compound)
all_sig$Group <- paste0(all_sig$Compound, "_", all_sig$Concentration)
all_sig$Passage <- as.integer(all_sig$Passage1)

p1_a <- ggplot(all_sig, aes(Passage, Fold_change, group = Group)) + 
  geom_point(aes(shape = FDR_sig), color = "gray74", size = 3, alpha = 0.6) +
  geom_point(data = all_sig[all_sig$Compound %in% sig_lst, ],
             aes(Passage, Fold_change, shape = FDR_sig, color = Compound),
             alpha = 0.9, size = 3) +
  scale_shape_manual(values = c(1, 16)) +
  geom_line(aes(linetype = Concentration), color = "gray", alpha = 0.6) + 
  geom_line(data = all_sig[all_sig$Compound %in% sig_lst, ],
             aes(Passage, Fold_change, 
                 group = Group, color = Compound, linetype = Concentration),
            alpha = 0.9) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  geom_hline(yintercept = 0, linetype = "twodash", 
             color = "gold", size = 0.8) +
  labs(y = "Fold change") +
  scale_color_manual(values = c("#39b87f", "#ba7e45", "#ea8783")) +
  #scale_color_manual(values = c("#94A323", "#EF8A0C", "#26897E")) +
  #geom_hline(yintercept = -0.2, linetype = "twodash", 
  #           color = "gold", size = 0.8) +
  main_theme #+
  #theme(legend.position = "top")

##############################################################################
ale1 <- read.table("../03.results/Figure_S3_growth_curve_24h_OD.txt",
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


p1_b <- ggplot(ale12_main[ale12_main$Compound == "DMSO_control", ],
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
  scale_color_manual(values = c("DMSO_control" = "gray74",
                                 "PFNA" = "lightblue4",
                                 "PFOA" = "yellow4")) 

p1 <- plot_grid(p1_a, p1_b, nrow = 2, align = "v", axis = "lr",
                rel_heights = c(1.2, 1), labels = c("a", "b"))

ggsave("../04.figures/Figure_1.pdf", p1, width = 10, height = 6.4)
