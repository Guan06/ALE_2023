### Make sure to run Figure_2.R before this scripts as that generated the filtered
### vcf files used here.

source("settings.R")
vcf <- readRDS("../04.results/20250910_evolved_population_variants_final.rds")

tab <- vcf %>% group_by(Compound, Concentration, Plate, Sample_ID) %>%
  summarise(Number_of_variants = n())

tab$Plate <- paste0("Plate_", tab$Plate)
tab$Concentration <- paste0(tab$Concentration, " µM")

p_s4a <- ggplot(tab, aes(Compound,
                         Number_of_variants)) + 
  #geom_boxplot(color = "gray47", outlier.shape = NA) + 
  geom_point(aes(color = Concentration),
            position = position_jitter(w = 0.1, h = 0), shape = 1) + 
  facet_grid(cols = vars(Plate),
             scales = "free", space = "free_x") +
  scale_color_manual(values = c("0 µM" = "gray",
                                "25 µM" = "#a0cbe8", 
                                "50 µM" = "navyblue"))  +
  scale_shape_manual(values = c(1, 2, 3, 5, 4)) + 
  main_theme +
  #ylim(0, 12) +
  labs(x = "", y = "Number of mutations") + 
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p_s4a
################################################################################

p_s4b  <- ggplot(tab, aes(Number_of_variants)) + 
  geom_bar(fill = "gray") +  theme_bw() + 
  labs(x = "Number of mutations", y = "Number of samples")
p_s4b


################################################################################


################################################################################

t <- tab %>% group_by(Compound, Concentration) %>% 
  summarise(Average = mean(Number_of_variants))

p_s4c <- ggplot(t, aes(Concentration, Average)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(color = "gray47", shape = 1) +
  geom_line(aes(group = Compound), linetype = "dashed", color = "gray47") +
  labs(y = "Mean number of mutations") +
  main_theme + theme_bw()

p_s4c
################################################################################

p_s4_bc <- plot_grid(p_s4b, p_s4c, nrow = 1, rel_widths = c(1, 1),
                     labels = c('b', 'c'))
p_s4 <- plot_grid(p_s4a, p_s4_bc, nrow = 2, rel_heights = c(1.5, 1),
                  labels = c('a', ''), align = "v", axis = "r")

ggsave("../05.figures/Figure_EV4.pdf", p_s4, width = 7, height = 7)

################################################################################

ggsave("../05.figures/revision_figure_2.pdf", p_s4c, width = 4, height = 3)