############################################################################
################## Get AUC via growthcurve 
############### Test FC and Sig of growth between treatment and controls
library(growthcurver)
dat_folder <- "../00.data/growth_curve_24h/"
files <- list.files(dat_folder, full.names = T, pattern = "tidy")
map <- read.table("../00.data/growth_curve_24h/date_to_passage_map.txt",
                  header = T)
well <- read.table("../00.data/growth_plate_layout_20220529.txt",
                   header = T, sep = "\t")
well_order <- c("Time", well$Well[1:96])
well$ID <- paste0(well$Well, "_Plate_", well$Plate)
well <- well[, c("ID", "Compound", "Concentration", "Solvent")]

all_dat <- c()

for (f in files){
  this_dat <- read.table(f, header = T, sep = "\t")
  
  this_date <- unlist(strsplit(basename(f), "_"))[2]
  this_passage <- map[map$Date == this_date, 2]
  plate_lst <- unique(this_dat$Plate)
  
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
  all_gc$Passage <- this_passage
  all_dat <- rbind(all_dat, all_gc)
}

all_dat$Passage <- as.character(all_dat$Passage)
all_dat$Plate <- as.character(all_dat$Plate)

all_dat$ID <- paste0(all_dat$sample, "_", all_dat$Plate)
all_dat <- merge(all_dat, well)

#### Calculate the growth AUC for ALE 2.0, test FC and sig 
dat2 <- read.table("../00.data/growth_curve_24h/ALE2_20230125_P0_P5_P10_P15_P20.txt",
                   header = T, sep = "\t")
well2<- read.table("../00.data/growth_plate_layout_20230125.txt",
                   header = T, sep = "\t")

### give Plate_5 to differ from ALE 1.0
well2$Plate <- "Plate_5"
well2_order <- c("Date", "Passage", "Time", well2$Well)

well2$ID <- paste0(well2$Well, "_", well2$Plate)
well2 <- well2[, c("ID", "Compound", "Concentration", "Solvent")]

dat2 <- dat2[, match(well2_order, colnames(dat2))]

all_gc2 <- c()
for (p in unique(dat2$Passage)) {
  this <- dat2[dat2$Passage == p, ]
  this <- this[, -c(1, 2)]
  
  this_gc <- SummarizeGrowthByPlate(this, plot_fit = F)
  this_gc$Passage <- p
  all_gc2 <- rbind(all_gc2, this_gc)
}
all_gc2$Plate <- "Plate_5"
all_gc2$ID <- paste0(all_gc2$sample, "_Plate_5")
all_gc2 <- merge(all_gc2, well2)

############################ Add ALE 2.0 data
all_dat <- rbind(all_dat, all_gc2)
write.table(all_dat, "../04.results/Figure_1_cd_growth_curve.txt",
            quote = F, sep = "\t", row.names = F)
