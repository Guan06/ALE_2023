## This script was used to generate figure of Passage 17 OD for reviewer #1's comments
source("plot_settings.R")
source("function_get_sig.R")

### Compare the growth curve
all_od <- read.table("../04.results/Figure_S3_growth_curve_24h_OD.txt",
                     header = T, sep = "\t")
p17_od <- all_od[all_od$Passage == "17" & all_od$Compound == "DMSO_control", ]

p0 <- ggplot(p17_od, aes(Time, OD, group = ID)) +
  geom_line(color = "gray42") +
  geom_point(shape = 1, size = 1.2, color = "gray42") +
  facet_wrap(~Plate, scales = "fixed", nrow = 1) +
  labs(x = "Time / h", y = "OD (600 nm)") +
  theme_bw()
p0

### Compare AUC of logistic curve
all_dat <- read.table("../04.results/Figure_1_cde_growth_curve.txt",
                      header = T, sep = "\t")

p17 <- all_dat[(all_dat$Passage == "17" & all_dat$Compound == "DMSO_control"), ]

p1 <- ggplot(p17, aes(Plate, auc_l)) +
  geom_boxplot(color = "gray42") +
  geom_jitter(shape = 1, color = "gray42") +
  labs(x = "", y = "AUC of logistic curve") +
  theme_bw()
p1

library(agricolae)
kruskal(p17$auc_l, p17$Plate, p.adj = "fdr")$groups

### Compare final OD of passage 17 in deep well plates
od <- read.table("../00.data/growth_OD595_20220529.txt", header = T, sep = "\t",
                 check.names = F)
well <- read.table("../00.data/growth_plate_layout_20220529.txt",
                   header = T, sep = "\t")

well$ID <- paste0(well$Well, "_", well$Plate)
well <- well[, c("ID", "Compound", "Concentration", "Solvent")]

od_long <- od %>% pivot_longer(
  !c("Passage", "Plate", "Row"),
  names_to = "Column",
  values_to = "OD")

od_long <- as.data.frame(od_long)
od_long$Well <- paste0(od_long$Row, od_long$Column)
od_long$ID <- paste0(od_long$Well, "_", od_long$Plate)

all_od <- merge(od_long, well)
all_od$Plate <- paste0("Plate_", all_od$Plate)

p17_od_595 <- all_od[all_od$Passage == "17" &
                       all_od$Compound == "DMSO_control", ]

p2 <- ggplot(p17_od_595, aes(Plate, OD)) +
  geom_boxplot(outlier.shape = NA, color = "gray42") +
  geom_jitter(shape = 1,color = "gray42") +
  labs(x = "", y = "OD (595 nm)") +
  theme_bw()
p2
kruskal(p17_od_595$OD, p17_od_595$Plate, p.adj = "fdr")$groups

p12 <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("b", "c"))

p <- plot_grid(p0, p12, nrow = 2, align = "v", axis = "r", labels = c("a", ""))

p
ggsave("../05.figures/revision_figure_1.pdf", p, width = 7, height = 4.5)
