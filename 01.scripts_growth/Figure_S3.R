source("plot_settings.R")

############################################################################
well <- read.table("../00.data/growth_plate_layout_20220529.txt",
                   header = T, sep = "\t")

well$ID <- paste0(well$Well, "_", well$Plate)
well <- well[, c("ID", "Compound", "Concentration", "Solvent")]

############################################################################
###################################### Panel a - Final OD between DMSO and 
###################################### Water control samples over 20 passages
od <- read.table("../00.data/growth_OD595_20220529.txt", header = T, sep = "\t",
                 check.names = F)
od_long <- od %>% pivot_longer(
            !c("Passage", "Plate", "Row"),
            names_to = "Column",
            values_to = "OD")

od_long <- as.data.frame(od_long)
od_long$Well <- paste0(od_long$Row, od_long$Column)
od_long$ID <- paste0(od_long$Well, "_", od_long$Plate)

all_od <- merge(od_long, well)
all_od$Plate <- as.character(all_od$Plate)
all_od$Passage <- as.character(all_od$Passage)

## correct the OD by subtracting the OD of non-contaminated mGAM
mgam <- all_od[all_od$Compound == "Negative_control", ]
mgam <- mgam[mgam$OD < 0.2, ]

mgam_mean <- mgam %>% group_by(Passage) %>%
    summarise(Mean_control = mean(OD))

all_od <- merge(all_od, mgam_mean)
all_od$OD_adjust <- all_od$OD - all_od$Mean_control

############################################################################
## Plot the OD of water and DMSO control over passage
controls <- all_od[all_od$Compound %in% c("Water_control", "DMSO_control"), ]
controls$Group <- paste0(controls$Passage, "_", controls$Compound)

## normalized the DMSO control by the mean OD of water control samples.
Water_control_mean <- controls %>% filter(Compound == "Water_control") %>%
  group_by(Passage) %>% summarise(Water_control_mean = mean(OD_adjust))

controls <- merge(controls, Water_control_mean)
controls$FC <- controls$OD_adjust / controls$Water_control_mean 

num <- length(unique(controls$Passage))
order <- paste0(rep(0:(num - 1), each = 2), "_",
                rep(c("Water_control", "DMSO_control"), num))

p_s3a <- ggplot(controls, aes(Group, FC, color = Compound)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(shape = Plate), alpha = 0.7, width = 0.05) +
    scale_x_discrete(limits = order) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    scale_shape_manual(values = c(1, 2, 3, 5)) +
    main_theme +
    labs(y = "Fold change") +
    geom_hline(yintercept = 1, color = "salmon", linetype = 2) +
    theme(axis.text.x = element_text(colour = "black", angle = 90,
              size = 6, hjust = 1))

############################################################################
## significant test
source("function_get_sig.R")
control_sig <- get_sig_control(controls, "Group", "OD_adjust")

write.table(control_sig, "../04.results/Figure_S3a_sig.txt", quote = F,
            row.names = F, sep = "\t")

############################################################################
###################################### Panel b, compare the growth curve 
###################################### of DMSO and Water control
dat_folder <- "../00.data/growth_curve_24h/"
files <- list.files(dat_folder, full.names = T, pattern = "tidy")
map <- read.table("../00.data/growth_curve_24h/date_to_passage_map.txt",
                  header = T)
all_dat <- c()

for (f in files){
  this_dat <- read.table(f, header = T, sep = "\t")
  
  this_date <- unlist(strsplit(basename(f), "_"))[2]
  this_passage <- map[map$Date == this_date, 2]
  
  this_dat_long <- this_dat %>% pivot_longer(!c(Plate, Time),
                                             names_to = "Well",
                                             values_to = "OD")
  this_dat_long$Passage <- this_passage
  all_dat <- rbind(all_dat, this_dat_long)
}
all_dat$ID <- paste0(all_dat$Well, "_", all_dat$Plate)

well <- read.table("../00.data/growth_plate_layout_20220529.txt",
                   header = T, sep = "\t")
well$ID <- paste0(well$Well, "_Plate_", well$Plate)
well <- well[, c("ID", "Compound", "Concentration", "Solvent")]

all_dat_meta <- merge(all_dat, well)
write.table(all_dat_meta, "../04.results/Figure_S3_growth_curve_24h_OD.txt",
            quote = F, sep = "\t", row.names = F)

############################################################################
## Plot the OD of water and DMSO control over passage
control_24h <- all_dat_meta[all_dat_meta$Compound %in% c("Water_control",
                                                         "DMSO_control"), ]

### filter out the outlier: H12_Plate_1 at Passage 0 at Time point 0
control_24h <- control_24h[!((control_24h$ID == "H12_Plate_1") &
                             (control_24h$Passage == 0) &
                               (control_24h$Time == 0)), ]
### filter out data <24h
control_24h <- control_24h[control_24h$Time < 24, ]

p_s3b <- ggplot(control_24h, aes(Time, OD, group = ID)) +
  geom_line(aes(color = Compound), alpha = 0.7) +
  geom_point(aes(color = Compound, shape = Plate),
             alpha = 0.5, size = 0.6) +
  facet_wrap(~Passage,
             scales = "free", nrow = 2) +
  main_theme +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_shape_manual(values = c(1, 2, 3, 5)) 

p_s3b_legend <- get_legend(p_s3b)

p_s3a <- p_s3a + theme(legend.position = "none")
p_s3b <- p_s3b + theme(legend.position = "none")

p_s3 <- plot_grid(p_s3a, p_s3b, nrow = 2, rel_heights = c(1, 1.5),
                  labels = c('a', 'b'))

ggsave("../05.figures/Figure_S3.pdf", p_s3, width = 10, height = 8)
ggsave("../05.figures/Figure_S3_legend.pdf", p_s3b_legend,
       width = 2, height = 3)
