### Make sure to run/source pre_Figure_S4.R before this scripts
### AF > 0.05
source("plot_settings.R")
tab1 <- read.table("../04.results/Figure_S4_variant_number_AF_005.txt",
                  header = T, sep = "\t")
tab1 <- tab1[tab1$Exp != "ALE2", ]

avg_nt1 <- mean(tab1[tab1$Compound == "NT5002", ]$Number_of_variants)
print(paste("Average number of variant (AF > 0.05) of parental strain: ",
            avg_nt1))
tab1$Group <- "AF > 0.05"

### AF > 0.5
tab2 <- read.table("../04.results/Figure_S4_variant_number_AF_05.txt",
                  header = T, sep = "\t")
tab2 <- tab2[tab2$Exp != "ALE2", ]
avg_nt2 <- mean(tab2[tab2$Compound == "NT5002", ]$Number_of_variants)
print(paste("Average number of variant (AF > 0.5) of parental strain: ",
            avg_nt2))
tab2$Group <- "AF > 0.5"

tab <- rbind(tab1, tab2)
tab$Concentration <- as.character(tab$Concentration)
tab$Plate <- as.character(tab$Plate)

df_avg <- data.frame(Group = unique((tab$Group)),
                     Avg = c(avg_nt1, avg_nt2))

p_s4a <- ggplot(tab, aes(Compound,
                         Number_of_variants)) + 
  geom_boxplot(color = "gray47", outlier.shape = NA) + 
  geom_point(aes(color = Concentration),
             alpha = 0.8,
             position = "jitter") + 
  facet_grid(rows = vars(Group), cols = vars(Plate),
             scales = "free", space = "free_x") +
  scale_color_manual(values = c("0" = "gray",
                                "25" = "#a0cbe8", 
                                "250" = "navyblue",
                                "50" = "#1170aa", 
                                "500" = "navyblue"))  +
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
  geom_bar(fill = "gray") +  main_theme + 
  geom_vline(xintercept = median_nt, linetype = "longdash", color = "salmon") +
  labs(x = "Number of variants (AF > 0.05)", y = "Number of samples")

print(paste("Median of non-parental", median_nt))

################################################################################
median_nt <- median(tab2[tab2$Compound != "NT5002", ]$Number_of_variants)
p_s4c  <- ggplot(tab2, aes(Number_of_variants)) + 
  geom_bar(fill = "gray") +  main_theme + 
  geom_vline(xintercept = median_nt, linetype = "longdash", color = "salmon") +
  labs(x = "Number of variants (AF > 0.5)", y = "Number of samples")

print(paste("Median of non-parental", median_nt))

################################################################################
p_s4_bc <- plot_grid(p_s4b, p_s4c, nrow = 1, labels = c('b', 'c'))

p_s4 <- plot_grid(p_s4a, p_s4_bc, nrow = 2, rel_heights = c(2, 0.8),
                  labels = c('a', ''), align = "v", axis = "r")

ggsave("../05.figures/Figure_S4.pdf", p_s4, width = 9, height = 8)
