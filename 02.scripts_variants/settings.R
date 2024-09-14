################################################################################
library(ggplot2)
library(cowplot)
library(stringr)
library(tidyr)
library(dplyr)

options(dplyr.summarise.inform = FALSE)

palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                      "#0072B2", "#D55E00", "#CC79A7", "#999999")

getOI = colorRampPalette(palette_OkabeIto)

main_theme <- theme(panel.background = element_blank(),
                    plot.background = element_blank(),
                    panel.grid = element_blank(),
                    axis.line.x = element_line(color = "black"),
                    axis.line.y = element_line(color = "black"),
                    axis.ticks = element_blank(),
                    axis.text = element_text(colour = "black", size = 10),
                    legend.background = element_blank())

################################################################################

concentration_color<-  c("0" = "gray",
                              "25" = "#a0cbe8", 
                              "250" = "navyblue",
                              "50" = "#1170aa", 
                              "500" = "navyblue")

effect_shape <- c("non-synonymous" = 16, "synonymous" = 1, "intergenic" = 3)


###############################################################################
## Colors and shapes in Figure 3
cc_color <- c("Parental0" = palette_OkabeIto[1],
              "Water0" = palette_OkabeIto[6],
              "Water_control0" = palette_OkabeIto[6],
              #"Xanthan_gum25" = "slateblue1",
              #"Xanthan_gum50" = "slateblue4")
              "Xanthan_gum25" = palette_OkabeIto[2],
              "Xanthan_gum50" = palette_OkabeIto[5])

cc_shape <- c("Parental" = 1,
              "Plate4_A4" = 1, "Plate4_A5" = 1, "Plate4_E4" = 1,
              "Plate4_B4" = 2, "Plate4_B5" = 2, "Plate4_F4" = 2,
              "Plate4_C4" = 5, "Plate4_C5" = 5, "Plate4_G4" = 5,
              "Plate4_D4" = 8, "Plate4_D5" = 8, "Plate4_H4" = 8)
