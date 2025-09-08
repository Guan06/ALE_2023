source("../02.scripts_variants/settings.R")

od_iso <- read.table("../00.data//20240912_iso_OD_parental_and_XG.txt",
                     header = T, sep = "\t")
od_iso$Date <- as.character(od_iso$Date)

### Remove obvious outlier caused by experimental error? e.g. the last row 
### (DMSO) was not inoculated, therefore no growth there
od_iso<- od_iso[!(od_iso$Row == "H" & od_iso$OD < 0.4), ]
od_iso <- od_iso[!(od_iso$Row == "D" & od_iso$OD > 0.4 & od_iso$Antibiotic == "Metronidazole"), ]

p_b <- ggplot(od_iso, aes(log2(Concentration), OD)) +
  geom_point(aes(shape = Date, color = Compound_conc2), alpha = 0.5) +
  scale_shape_manual(values = c(1, 2, 3, 5, 8)) +
  #scale_shape_manual(values = c(1, 1, 1, 2, rep(3, 6), 4, rep(9, 7))) +
  geom_line(aes(group = Group, color = Compound_conc2), alpha = 0.3) +
  facet_grid(cols = vars(Antibiotic), rows = vars(Compound_conc2), 
             scales = "free") +
  #facet_wrap(~Antibiotic * Evol_compound, scales = "free_x", nrow = 2) +
  main_theme +
  #scale_x_reverse() +
  scale_color_manual(values = cc_color) +
  geom_hline(yintercept = 1, color = "gray", linetype = "dashed") +
  geom_hline(yintercept = 0.4, color = "gray", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "top")

ggsave("../05.figures/Figure_EV8.pdf", p_b, height = 6, width = 10)

################################################################################
############################# Facet each population 
od_iso <- od_iso %>% mutate(Rep = case_when(
  Date == "20240621" ~ 'Parental',
  (Date == "20240720" | Date == "20240705") ~ 'Batch 1',
  (Date == "20240802" | Date == "20240808") ~ 'Batch 2',
))

od_iso$Population <- factor(od_iso$Population,
                            levels = c("Parental",
                                       "Plate4_A5", "Plate4_B5", "Plate4_C5",
                                       "Plate4_A4", "Plate4_B4", "Plate4_C4", "Plate4_D4",
                                       "Plate4_E4", "Plate4_F4", "Plate4_G4", "Plate4_H4"))

p_ev8_v2 <- ggplot(od_iso, aes(log2(Concentration), OD)) +
  geom_point(aes(shape = Rep, color = Compound_conc2), alpha = 0.5) +
  geom_line(aes(group = Group, color = Compound_conc2), alpha = 0.3) +
  facet_grid(cols = vars(Antibiotic), rows = vars(Population), 
             scales = "free") +
  #main_theme +
  theme_minimal_vgrid() + 
  scale_shape_manual(values = c(2, 1, 5)) +
  scale_color_manual(values = cc_color) +
  geom_hline(yintercept = 1, color = "gray", linetype = "dashed") +
  geom_hline(yintercept = 0.4, color = "gray", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "top")

p_ev8_v2
ggsave("../05.figures/Figure_EV8_v2.pdf", p_ev8_v2, height = 12, width = 8)
