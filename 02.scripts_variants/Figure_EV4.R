### Make sure to run/source pre_Figure_S4.R before this scripts
### AF > 0.05
source("settings.R")
tab1 <- read.table("../04.results/Figure_EV4_variant_number_AF_005.txt",
                  header = T, sep = "\t")
tab1 <- tab1[tab1$Exp != "ALE2", ]

avg_nt1 <- mean(tab1[tab1$Compound == "NT5002", ]$Number_of_variants)
print(paste("Average number of variant (AF > 0.05) of parental strain: ",
            avg_nt1))
tab1$Group <- "AF > 0.05"

### AF > 0.5
tab2 <- read.table("../04.results/Figure_EV4_variant_number_AF_05.txt",
                  header = T, sep = "\t")
tab2 <- tab2[tab2$Exp != "ALE2", ]
avg_nt2 <- mean(tab2[tab2$Compound == "NT5002", ]$Number_of_variants)
print(paste("Average number of variant (AF > 0.5) of parental strain: ",
            avg_nt2))
tab2$Group <- "AF > 0.5"

tab <- rbind(tab1, tab2)
df_avg <- data.frame(Group = unique((tab$Group)),
                     Avg = c(avg_nt1, avg_nt2))

tab$Plate <- paste0("Plate_", tab$Plate)
tab$Concentration <- paste0(tab$Concentration, " µM")

p_s4a <- ggplot(tab, aes(Compound,
                         Number_of_variants)) + 
  geom_boxplot(color = "gray47", outlier.shape = NA) + 
  geom_point(aes(color = Concentration),
             position = "jitter", shape = 1) + 
  facet_grid(rows = vars(Group), cols = vars(Plate),
             scales = "free", space = "free_x") +
  scale_color_manual(values = c("0 µM" = "gray",
                                "25 µM" = "#a0cbe8", 
                                "50 µM" = "navyblue"))  +
  scale_shape_manual(values = c(1, 2, 3, 5, 4)) + 
  geom_abline(data = df_avg, aes(intercept = Avg, slope = 0),
              linetype = "longdash", color = "salmon") +
  main_theme +
  labs(x = "", y = "Number of variants") + 
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

################################################################################
median_nt <- median(tab1[tab1$Compound != "NT5002", ]$Number_of_variants)
p_s4b  <- ggplot(tab1, aes(Number_of_variants)) + 
  geom_bar(fill = "gray") +  theme_bw() + 
  geom_vline(xintercept = median_nt, linetype = "longdash", color = "salmon") +
  labs(x = "Number of variants (AF > 0.05)", y = "Number of samples")

print(paste("Median of non-parental", median_nt))

################################################################################
median_nt <- median(tab2[tab2$Compound != "NT5002", ]$Number_of_variants)
p_s4c  <- ggplot(tab2, aes(Number_of_variants)) + 
  geom_bar(fill = "gray") +  theme_bw() + 
  geom_vline(xintercept = median_nt, linetype = "longdash", color = "salmon") +
  labs(x = "Number of variants (AF > 0.5)", y = "Number of samples")

print(paste("Median of non-parental", median_nt))

################################################################################

tab_evolved <- tab[tab$Compound != "NT5002", ]
t <- tab_evolved %>% group_by(Group, Compound, Concentration) %>% 
  summarise(Average = mean(Number_of_variants))

p_s4d <- ggplot(t, aes(Concentration, Average)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(color = "gray47", shape = 1) +
  geom_line(aes(group = Compound), linetype = "dashed", color = "gray47") +
  facet_wrap(~ Group, scales = "free_y") +
  labs(y = "Mean number of detected variants") +
  main_theme + theme_bw()
################################################################################

p_s4_bcd <- plot_grid(p_s4b, p_s4c, p_s4d, nrow = 1,
                      rel_widths = c(1, 1, 1.4), labels = c('b', 'c', 'd'))
p_s4 <- plot_grid(p_s4a, p_s4_bcd, nrow = 2, rel_heights = c(2, 1.1),
                  labels = c('a', ''), align = "v", axis = "r")

ggsave("../05.figures/Figure_EV4.pdf", p_s4, width = 9, height = 8)

################################################################################

ggsave("../05.figures/revision_figure_2.pdf", p_s4d, width = 5, height = 4)