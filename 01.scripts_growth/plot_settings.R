library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(tidyr)
library(dplyr)
library(ggpubr)

getDark2 = colorRampPalette(brewer.pal(8, "Dark2"))
getSet1 = colorRampPalette(brewer.pal(9, "Set1"))

BrBG <- brewer.pal(n = 11, name = "BrBG")

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
                    text = element_text(family="sans"))
