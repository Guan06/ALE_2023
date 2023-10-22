### Make sure to run/source pre_Figure_S5.R before this scripts

tab <- read.table("../03.results/Figure_S5_variant_number_AF_005.txt",
                  header = T, sep = "\t")
tab$Concentration <- as.character(tab$Concentration)
tab$Plate <- as.character(tab$Plate)
avg_nt <- mean(tab[tab$Compound == "NT5002", ]$Number_of_variants)
median_nt <- median(tab[tab$Compound != "NT5002", ]$Number_of_variants)

p_s5a <- ggplot(tab, aes(Compound,
                         Number_of_variants)) + 
  geom_boxplot(color = "gray47", outlier.shape = NA) + 
  geom_point(aes(color = Concentration, shape = Plate),
             alpha = 0.8,
             position = "jitter") + 
  facet_grid(~Plate, space = "free", scales = "free_x") +
  scale_color_manual(values = c("0" = "gray",
                                "25" = "#a0cbe8", 
                                "250" = "navyblue",
                                "50" = "#1170aa", 
                                "500" = "navyblue"))  +
  scale_shape_manual(values = c(1, 2, 3, 5, 4)) + 
  geom_hline(yintercept = avg_nt, linetype = "longdash", color = "salmon") +
  main_theme +
  labs(x = "", y = "No. of variants (AF > 0.05)") + 
  theme(legend.position = "top",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  

p_s5c  <- ggplot(tab, aes(Number_of_variants)) + 
  geom_bar(fill = "gray") +  main_theme + 
  geom_vline(xintercept = median_nt, linetype = "longdash", color = "salmon") +
  labs(x = "Number of variants (AF > 0.05)", y = "Number of samples")

print(paste("Mean of parental:", avg_nt, "Median of non-parental", median_nt))

################################################################################

tab <- read.table("../03.results/Figure_S5_variant_number_AF_05.txt",
                  header = T, sep = "\t")
tab$Concentration <- as.character(tab$Concentration)
tab$Plate <- as.character(tab$Plate)
avg_nt <- mean(tab[tab$Compound == "NT5002", ]$Number_of_variants)
median_nt <- median(tab[tab$Compound != "NT5002", ]$Number_of_variants)

p_s5b <- ggplot(tab, aes(Compound,
                         Number_of_variants)) + 
  geom_boxplot(color = "gray47", outlier.shape = NA) + 
  geom_point(aes(color = Concentration, shape = Plate),
             alpha = 0.8,
             position = "jitter") + 
  labs(x = "", y = "No. of variants (AF > 0.5)") + 
  facet_grid(~Plate, space = "free", scales = "free") +
  scale_color_manual(values = c("0" = "gray",
                                "25" = "#a0cbe8", 
                                "250" = "navyblue",
                                "50" = "#1170aa", 
                                "500" = "navyblue"))  +
  scale_shape_manual(values = c(1, 2, 3, 5, 4)) + 
  geom_hline(yintercept = avg_nt, linetype = "longdash", color = "salmon") +
  main_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

p_s5d  <- ggplot(tab, aes(Number_of_variants)) + 
  geom_bar(fill = "gray") +  main_theme + 
  geom_vline(xintercept = median_nt, linetype = "longdash", color = "salmon") +
  labs(x = "Number of variants (AF > 0.5)", y = "Number of samples")

print(paste("Mean of parental:", avg_nt, "Median of non-parental", median_nt))

################################################################################
p_s5_ab <- plot_grid(p_s5a, p_s5b, nrow = 2, rel_heights = c(1, 1.25),
                     labels = c('a', 'b'))
p_s5_cd <- plot_grid(p_s5c, p_s5d, nrow = 1, labels = c('c', 'd'))

p_s5 <- plot_grid(p_s5_ab, p_s5_cd, nrow = 2, rel_heights = c(2, 0.8))
ggsave("../04.figures/Figure_S5.pdf", p_s5, width = 9, height = 8)
