source("settings.R")

od_pop <- read.table("../00.data/20240912_pop_OD_parental_and_XG.txt",
                     header = T, sep = "\t")
colnames(od_pop)[10] <- "Compound"
od_pop$Group2 <- paste0(od_pop$Compound, od_pop$Compound_conc)

ggplot(od_pop, aes((Concentration), OD)) +
  geom_point(aes(shape = Rep, color = Group2), alpha = 0.6) +
  scale_shape_manual(values = c(1, 2, 3)) +
  #scale_shape_manual(values = c(1, 1, 1, 2, rep(3, 6), 4, rep(9, 7))) +
  geom_line(aes(group = Group, color = Group2), alpha = 0.5) +
  facet_grid(cols = vars(Antibiotic), rows = vars(Compound), 
             scales = "free") +
  #facet_wrap(~Antibiotic * Evol_compound, scales = "free_x", nrow = 2) +
  main_theme +
  #scale_x_reverse() +
  scale_color_manual(values = cc_color) +
  geom_hline(yintercept = 1, color = "gray", linetype = "dashed") +
  geom_hline(yintercept = 0.4, color = "gray", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "top")

ggsave("../04.figures/Figure_S8.pdf", height = 6, width = 10)
