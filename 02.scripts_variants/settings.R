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
                    axis.text = element_text(colour = "black", size = 10))
################################################################################

concentration_color<-  c("0" = "gray",
                              "25" = "#a0cbe8", 
                              "250" = "navyblue",
                              "50" = "#1170aa", 
                              "500" = "navyblue")

c
effect_shape <- c("non-synonymous" = 16, "synonymous" = 1, "intergenic" = 3)


