source("settings.R")

###############################################################################
# read in all variants detected in ALE 1.0 and ALE 2.0 with AF > 0.5
all <- readRDS("../03.results/20231023_all_05.rds")
t1 <- all %>% group_by(POS) %>% summarise(Prevalence =
                                            length(unique(Sample_ID)),
                                          Average_AF = mean(Ratio))
t1_des <- unique(all[, c("POS", "FTYPE", "Effect")])
t1 <- merge(t1, t1_des)

map <- read.table("../00.data/vcf_20231024_effect_map.txt", 
                  header = T, sep = "\t")
t1 <- merge(t1, map)
p2_1 <- ggplot(t1, aes(POS, Prevalence)) +
  geom_point(aes(color = FTYPE, size = Average_AF, shape = Effect_type),
             alpha = 0.7) + 
  scale_color_manual(values = c("gold", "gray47")) +
  xlab("Position") +
  scale_shape_manual(values = map$Shape) +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  scale_size(range = c(2.5, 5)) +
  xlim(c(0, 5e+06)) +
  guides(colour = guide_legend(order = 1, nrow = 2), 
         size = guide_legend(order = 2, nrow = 2),
         shape = guide_legend(order = 3, nrow = 2)) +
  theme(legend.position = "bottom",
        legend.justification = c("left", "center")) +
  main_theme 

p2_2 <- ggplot(all, aes(POS)) + 
  geom_density(aes(color = FTYPE, y = ..scaled..)) +
  scale_color_manual(values = c("gold", "gray47")) +
  scale_y_continuous(position = "right") +
  xlab("") + ylab("Density of variants") +
  xlim(c(0, 5e+06)) +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  theme(legend.position = "none",
        axis.text.x=element_blank()) +
  main_theme

aligned_plots <- align_plots(p2_2, p2_1, align="hv", axis="tblr")
p2 <- ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])

###############################################################################
#p <- plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 1.25), 
#               align = "v", axis = "lr", labels = c("a", "b"))
ggsave("../04.figures/Figure_S6.pdf", p2, width = 8, height = 2.5)
