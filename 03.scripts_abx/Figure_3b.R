source("settings.R")
library(ggridges)

od_iso <- read.table("../00.data//20240912_iso_OD_parental_and_XG.txt",
                     header = T, sep = "\t")
od_iso$Date <- as.character(od_iso$Date)

### Remove obvious outlier caused by experimental error? e.g. the last row 
### (DMSO) was not inoculated, therefore no growth there
od_iso<- od_iso[!(od_iso$Row == "H" & od_iso$OD < 0.4), ]
od_iso <- od_iso[!(od_iso$Row == "D" & od_iso$OD > 0.4 & od_iso$Antibiotic == "Metronidazole"), ]

t1 <- od_iso %>% group_by(Group, Antibiotic) %>% 
  mutate(Min_OD = min(OD), Max_OD = max(OD)) %>% 
  mutate(MIC80_OD = (Max_OD - Min_OD) * 0.2 + Min_OD)


t1 <- t1[t1$OD > t1$MIC80_OD, ]

t1 <- t1 %>% group_by(Group, Antibiotic) %>% mutate(MIC80 = max(Concentration))

t2 <- unique(t1[, c("Group", "Compound", "Compound_conc2", 
                    "MIC80", "Antibiotic", "Date", "Population")])

t2 <- t2 %>% mutate(Rep = case_when(
  Date == "20240621" ~ 'Parental',
  (Date == "20240720" | Date == "20240705") ~ 'Rep_1',
  (Date == "20240802" | Date == "20240808") ~ 'Rep_2',
))

ggplot(t2, aes(y = interaction(Compound_conc2, Rep), x = MIC80,
               fill = Compound_conc2)) +
  geom_jitter(aes(color = Compound_conc2, shape = Population), alpha = 0.6,
              width = 0, height = 0.3) +
  geom_density_ridges(scale = 0.5, quantile_lines = F, quantiles = 2,
                      color = "gray67", lwd = 0.3, alpha = 0.7) +
  facet_grid(cols=vars(Antibiotic), rows=vars(Rep),
             scales = "free", space = "free_y") +
  scale_color_manual("Population (compound and concentration)", values = cc_color) +
  scale_fill_manual(values = cc_color, guide = F) +
  scale_shape_manual(values = cc_shape) +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape = guide_legend(order = 2, nrow = 2)) +
  #  scale_shape_manual(values = c(5, 2, 1)) +
  main_theme +
  theme(axis.text.y = element_blank(),
        legend.position = "top") +
  labs(y = "Compound and concentration", x = "MIC")

ggsave("../05.figures/Figure_3b.pdf", width = 10, height = 6)
