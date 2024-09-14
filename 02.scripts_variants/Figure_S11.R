source("settings.R")

od_iso <- read.table("../00.data//20240912_iso_OD_parental_and_XG.txt",
                     header = T, sep = "\t")
od_iso$Date <- as.character(od_iso$Date)

p_b <- ggplot(od_iso, aes(log(Concentration + 0.01), OD)) +
  geom_point(aes(shape = Date, color = Compound_conc2), alpha = 0.5) +
  scale_shape_manual(values = c(1, 2, 3, 5, 8)) +
  #scale_shape_manual(values = c(1, 1, 1, 2, rep(3, 6), 4, rep(9, 7))) +
  geom_line(aes(group = Group, color = Compound_conc2), alpha = 0.3) +
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

p <- plot_grid(p_a, p_b, nrow = 1, labels = c("a", "b"), rel_widths = c(1, 3))

ggsave("../04.figures/Figure_S11.pdf", p, height = 6, width = 10)
