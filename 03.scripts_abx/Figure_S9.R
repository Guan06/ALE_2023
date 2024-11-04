source("settings.R")

#############################################  parental strain
od_des <- read.table("../00.data/20240912_iso_OD_parental_and_XG.txt",
                     header = T, sep = "\t")
od_des3 <- od_des[od_des$Population == "Parental", ]
stat <- od_des3 %>% group_by(Row, Antibiotic, Concentration) %>% 
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

data_df <- stat[, c("Row_new", "Column", "Mean_OD", "Concentration")]

# Create the plot
p <- ggplot(data_df, aes(x = Column, y = Row_new)) +
  geom_point(aes(fill = Mean_OD), size = 20, shape = 21,
             color = "gray") +  
  geom_text(aes(label = round(Concentration, 2)), color = "gray28",
            vjust = 0.5, size = 4.6) + 
  scale_fill_gradient("Avg.OD", low = "deepskyblue3", high = "orange", 
                      limits = c(0.15, 1.05)) +
  scale_size_identity()  +  # Ensure circles are of uniform size
  theme_void() +  # Remove axes and background
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +  # Add padding to x-axis
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.1))) +  # Add padding to y-axis
  theme(
    panel.grid = element_blank(),  
    legend.position = "right" 
  )

# Add column names (1 to 6)
p <- p + geom_text(data = data.frame(Column = 1:6,
                                     Row = 8.5), 
                   aes(x = Column, y = Row, label = c("Amoxicillin",
                                                      "Ampicillin",
                                                      "Chloramphenicol",
                                                      "Doxycycline",
                                                      "Erythromycin",
                                                      "Metronidazole")), 
                   angle = 30, vjust = 0, hjust = 0, size = 5,
                   inherit.aes = FALSE)

# Add row names (A to H)
p <- p + geom_text(data = data.frame(Row = 1:8, Column = 0.5), 
                   aes(x = Column, y = 9 - Row, label = LETTERS[1:8]), 
                   hjust = 1, size = 5, inherit.aes = FALSE)

ggsave("../05.figures/Figure_S9.pdf", p, width = 6, height = 6)
