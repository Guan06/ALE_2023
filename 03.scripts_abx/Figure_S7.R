source("settings.R")

od_pop <- read.table("../00.data/20241120_pop_OD_all_Abx_parental_and_XG.txt",
                     header = T, sep = "\t")
colnames(od_pop)[10] <- "Compound"
od_pop$Group2 <- paste0(od_pop$Compound, od_pop$Compound_conc)

##Filter out negative control
od_pop <- od_pop[od_pop$Antibiotic != "Negative_control", ]

ggplot(od_pop, aes(log2(Concentration), OD)) +
  geom_point(aes(shape = Rep, color = Group2), alpha = 0.6) +
  scale_shape_manual(values = c(1, 2, 3)) +
  geom_line(aes(group = Group, color = Group2), alpha = 0.5) +
  facet_grid(cols = vars(Antibiotic), rows = vars(Compound), 
             scales = "free") +
  main_theme +
  #scale_x_reverse() +
  scale_color_manual(values = cc_color) +
  geom_hline(yintercept = 1, color = "gray", linetype = "dashed") +
  geom_hline(yintercept = 0.4, color = "gray", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "top")

ggsave("../05.figures/Figure_S7.pdf", height = 6, width = 12)
