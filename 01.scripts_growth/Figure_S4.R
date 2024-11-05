source("plot_settings.R")
all_dat <- read.table("../04.results/Figure_1_cd_growth_curve.txt",
                      header = T, sep = "\t")

lma_lst <- c("Ezetimibe", "Rosuvastatin", "Simvastatin", "DMSO_control")
lma_dat <- all_dat[all_dat$Compound %in% lma_lst, ]
## only take DMSO_control from Plate_2 and Plate_3
lma_dat <- lma_dat[lma_dat$Plate %in% c("Plate_2", "Plate_3"), ]

lma_dat$Passage <- as.integer(lma_dat$Passage)
lma_dat$Concentration <- as.character(lma_dat$Concentration)

p_s4 <- ggplot(lma_dat, aes(interaction(Concentration,Compound, Plate, Passage),
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
  scale_alpha_manual(values = c("0" = 0.7,
                                "25" = 0.1,
                                "50" = 0.6)) +
  scale_color_manual(values = c("DMSO_control" = "gray",
                                "Ezetimibe"= "#39b87f",
                                "Rosuvastatin" = "#b07aa1", 
                                "Simvastatin"= "#ea8783")) +
  scale_fill_manual(values = c("DMSO_control" = "gray99",
                                "Ezetimibe"= "#39b87f",
                                "Rosuvastatin" = "#b07aa1", 
                                "Simvastatin"= "#ea8783"), guide = "none") +
  scale_shape_manual(values = c(2, 3)) +
  labs(y = "AUC of logistic curve", 
       x = "Passage 0 -> 4 -> 7 -> 10 -> 14 -> 17 -> 20") +
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 3),
         alpha = guide_legend(order = 2)) +
  theme(legend.position = "top", axis.text.x = element_blank())

ggsave("../05.figures/Figure_S4.pdf", p_s4, width = 8, height = 3)

# lma_dat$Group <- paste0(lma_dat$Concentration, "_", lma_dat$Compound,
#                         "_", lma_dat$Plate, "_", lma_dat$Passage)
# source("function_get_sig.R")
# lma_cmp <- get_sig(lma_dat, "Group", auc_l)

all_sig <- read.table("../04.results/Figure_1c_all_gc_cmp.txt", header = T,
                   sep = "\t")
lma_sig <- all_sig[all_sig$Compound %in% lma_lst, ]
write.table(lma_sig, "../04.results/Figure_S4_sig.txt",
            quote = F, row.names = F, sep = "\t")
