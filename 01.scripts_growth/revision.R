
## Scripts to process and analyse new data from the experiment during revision
library(colorblindcheck)
palette_check(c("#98dd94", "#E69F00", "#808060"), plot = TRUE)

library(growthcurver)
## Plate 1 and Plate 2 contains mGAM + DMSO
## Plate 3 and Plate 4 contains mGAM + DMSO + 50 µM Loperamide

this_dat <- read.table("../00.data/growth_curve_24h/tidy_20250514_Buni_Loperamide.txt",
                      header = T, sep = "\t")
well <- read.table("../00.data/growth_plate_layout_20250515.txt", 
                     header = T, sep = "\t")

well_order <- c("Time", well$Well[1:96])

plate_lst <- unique(this_dat$Plate)
plate_df <- data.frame(Plate = c("Plate_1", "Plate_2", "Plate_3", "Plate_4"),
                       Condition = c("DMSO", "DMSO",
                                     "Loperamide", "Loperamide"))

all_gc <- c()
for (pl in plate_lst){
  this_pl <- this_dat[this_dat$Plate == pl, ]
  this_pl <- this_pl[, -1]
  this_pl <- this_pl[, match(well_order, colnames(this_pl))]
  
  if (nrow(this_pl) >24) {this_pl <- this_pl[1:24, ]}
  
  this_gc <- SummarizeGrowthByPlate(this_pl, plot_fit = F)
  this_gc$Plate <- pl
  all_gc <- rbind(all_gc, this_gc)
}

colnames(all_gc)[1] <- "Well"
all_dat <- merge(all_gc, well)
all_dat_noNC <- all_dat[all_dat$Isolate != "NC", ] 
all_dat_noNC <- merge(all_dat_noNC, plate_df)

### Filter out Plate 1 G6 well as we saw some unsolvable material in that well
### Could be either paper or plastic pieces
all_dat_noNC <- all_dat_noNC[!(all_dat_noNC$Plate == "Plate_1" &
                                all_dat_noNC$Well == "G6"), ]
all_dat_noNC$Row_num <- as.factor(all_dat_noNC$Row_num)
all_dat_noNC$Group <- factor(all_dat_noNC$Group,
                                levels = c("Parental", "DMSO_evo", "Lop50_evo"),
                                ordered = TRUE)
p1e_right <- ggplot(all_dat_noNC, 
                    aes(interaction(Row_num, Group), auc_l)) +
  # geom_boxplot(aes(color = Group),
  #              outlier.shape = NA,
  #              position = position_dodge2(preserve = "single")) + 
  geom_boxplot(aes(color = Group), outlier.shape = NA) + 
  geom_point(aes(color = Group),#, shape = Row_num),
             position = position_jitterdodge(), alpha = 1, shape = 1) +
  facet_wrap(~Condition) +
  #scale_x_discrete(limits = c("Parental", "DMSO_evo", "Lop50_evo")) +
  scale_color_manual(values = c("Lop50_evo" = "#98dd94",
                                "Parental" = "#808060",
                                "DMSO_evo" = "#E69F00")) +
  #scale_shape_manual(values = 1:8) +
  main_theme +
  theme(axis.text.x = element_text(colour = "black", angle = 90,
                               size = 6, hjust = 1),
        legend.position = "NA") +
  theme(axis.text.x = element_blank()) +
  labs(y = "AUC of logistic curve", 
       x = "Isolates 1 to 8 of parental strain, DMSO- and Loperamide-evolved strain") 

p1e_right

##### Plot the growth curve
od <- this_dat %>% pivot_longer(!c("Plate", "Time"), 
                                names_to = "Well", values_to = "OD")
od <- merge(od, well)
od <- od[od$Group != "NC", ]

### Filter out Plate 1 G6 well as we saw some unsolvable material in that well
### Could be either paper or plastic pieces
od <- od[!(od$Plate == "Plate_1" & od$Well == "G6"), ]
od <- merge(od, plate_df)
od$ID2 <- paste0(od$Plate, "_", od$ID) 

p1e_left <- ggplot(od, aes(Time, OD, group = ID2)) +
  geom_line(aes(color = Group), alpha = 0.7) +
  geom_point(aes(color = Group, shape = Plate),
             alpha = 0.5, size = 1.2) +
  facet_wrap(~Condition, scales = "fixed", nrow = 1) +
  main_theme +
  scale_color_manual(values = c("Lop50_evo" = "#98dd94",
                                "Parental" = "#808060",
                                "DMSO_evo" = "#E69F00")) +
  scale_shape_manual(values = c(1, 2, 3, 5))  

p1e_left

p1e <- plot_grid(p1e_left, p1e_right, nrow = 1, labels = c("e", ""),
                align = "h", rel_widths = c(2, 1))
p1e

ggsave("../05.figures/revision.pdf", width = 12, height = 5)
