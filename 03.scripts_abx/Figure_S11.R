source("settings.R")

od_iso <- read.table("../00.data/20240912_iso_OD_parental_and_XG.txt",
                     header = T, sep = "\t")
od_iso$Date <- as.character(od_iso$Date)

#################################################################################
### Remove obvious outlier caused by experimental error? e.g. the last row 
### (DMSO) was not inoculated, therefore no growth there
od_iso<- od_iso[!(od_iso$Row == "H" & od_iso$OD < 0.4), ]
od_iso <- od_iso[!(od_iso$Row == "D" & od_iso$OD > 0.4 & od_iso$Antibiotic == "Metronidazole"), ]

### Panel A, split circle of all avg. OD
od_iso_no_parental <- od_iso[od_iso$Compound != "Parental", ]

#od_iso_no_parental %>% group_by(Row, Antibiotic, Compound) %>%
#  summarise(n())

stat <- od_iso_no_parental %>% group_by(Row, Antibiotic, Compound) %>% 
  summarise(Mean_OD = mean(OD))

map_col <- data.frame(Column = 1:6,
                      Antibiotic = c("Amoxicillin",
                                     "Ampicillin",
                                     "Chloramphenicol",
                                     "Doxycycline",
                                     "Erythromycin",
                                     "Metronidazole"))
map_row <- data.frame(Row_new = factor(1:8, levels = 8:1),
                      Row = LETTERS[1:8])

stat <- merge(stat, map_col)
stat <- merge(stat, map_row)

######### Function definition - Function to create half circles
half_circle <- function(x, y, radius, start_angle, end_angle) {
  tibble(
    x = x + radius * cos(seq(start_angle, end_angle, length.out = 100)),
    y = y + radius * sin(seq(start_angle, end_angle, length.out = 100))
  )
}
#########

stat <- stat[, c("Row_new", "Column", "Mean_OD", "Compound")]
stat <- stat %>% 
  pivot_wider(names_from = "Compound", values_from = "Mean_OD")

# Generate data for left and right halves of each well
left_half_data <- do.call(rbind, mapply(function(x, y, fill) {
  circle_data <- half_circle(x, y, radius = 0.4, start_angle = pi/2, end_angle = 3*pi/2)
  cbind(circle_data, fill = fill, group = paste(x, y, "left"))
}, stat$Column, as.numeric(stat$Row_new), stat$Water_control, SIMPLIFY = FALSE))

right_half_data <- do.call(rbind, mapply(function(x, y, fill) {
  circle_data <- half_circle(x, y, radius = 0.4, start_angle = -pi/2, end_angle = pi/2)
  cbind(circle_data, fill = fill, group = paste(x, y, "right"))
}, stat$Column, as.numeric(stat$Row_new), stat$Xanthan_gum, SIMPLIFY = FALSE))

# Combine the data
plot_data <- rbind(left_half_data, right_half_data)
plot_data <- as.data.frame(plot_data)  # Ensure it's a data frame

# Base plot
plot_data$y <- plot_data$y - 1.5

p_a <- ggplot(plot_data, aes(x = x, y = y, fill = fill, group = group)) +
  geom_polygon(color = "gray47", alpha = 1, lwd = 0.8) +
  scale_fill_gradient("Avg.OD", low = "deepskyblue3", high = "orange", 
                      limits = c(0.15, 1.05)) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(legend.position = "right") +
  theme(
    panel.grid = element_blank()
  )

# Add column names (1 to 6)
p_a <- p_a + geom_text(data = data.frame(Column = 1:6,
                                         Row = 7),
                       aes(x = Column, y = Row, label = c("Amoxicillin",
                                                          "Ampicillin",
                                                          "Chloramphenicol",
                                                          "Doxycycline",
                                                          "Erythromycin",
                                                          "Metronidazole")), 
                       angle = 30, vjust = 0, hjust = 0, size = 5,
                       inherit.aes = FALSE)

# Add row names (A to H)
p_a <- p_a + geom_text(data = data.frame(Row = 1:8, Column = 0.2), 
                       aes(x = Column, y = 7.5 - Row, label = LETTERS[1:8]), 
                       hjust = 1, size = 5, inherit.aes = FALSE)

p_a <- p_a + expand_limits(x = c(0, 7.4), y = c(-1, 8))

ggsave("../05.figures/Figure_S11.pdf", p_a, width = 7, height = 7)

## Add significance test for each well
od_iso_no_parental$Compare <- paste0(od_iso_no_parental$Row, "_",
                                     od_iso_no_parental$Antibiotic)
all <- c()
for (cmp in unique(od_iso_no_parental$Compare)) {
  this_cmp <- od_iso_no_parental[od_iso_no_parental$Compare == cmp, ]
  #ggplot(this_cmp, aes(Compound, OD)) + geom_violin()
  
  this_cmp1 <- this_cmp[this_cmp$Compound == "Water_control", ]
  this_cmp2 <- this_cmp[this_cmp$Compound == "Xanthan_gum", ]
  
  median_1 <- median(this_cmp1$OD)
  median_2 <- median(this_cmp2$OD)
  this_p <- wilcox.test(this_cmp1$OD, this_cmp2$OD)$p.value
  
  all <- rbind(all, c(cmp, median_1, median_2, this_p))
}

colnames(all) <- c("Comparison", 
                   "Median_OD_Water_control", "Median_OD_Xanthan_gum",
                   "P_value")

all <- as.data.frame(all)
all$P_adjust <- p.adjust(all$P_value, method = "bonferroni")
all$Sig <- ifelse(all$P_adjust < 0.01, "Yes", "No")
write.table(all, "../04.results/Figure_S11_sig.txt", quote = F, sep = "\t",
            row.names = F)
