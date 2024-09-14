source("settings.R")

###############################################################################
# read in all variants detected in ALE 1.0 and ALE 2.0 with AF > 0.05
all <- readRDS("../03.results/20231023_all_005.rds")
t1 <- all %>% group_by(POS) %>% summarise(Prevalence =
                                            length(unique(Sample_ID)),
                                          Average_AF = mean(Ratio))
t1_des <- unique(all[, c("POS", "FTYPE", "Effect")])
t1 <- merge(t1, t1_des)

map <- read.table("../00.data/vcf_20231024_effect_map.txt", 
                  header = T, sep = "\t")

t1 <- merge(t1, map)

p1_1 <- ggplot(t1, aes(POS, Prevalence)) +
  geom_point(aes(color = FTYPE, size = Average_AF, shape = Effect_type),
             alpha = 0.7) + 
  scale_shape_manual(limits = map$Effect_type,
                     values = map$Shape) +
  scale_color_manual(values = c("gold", "gray47")) +
  xlim(c(0, 5e+06)) + xlab("") +
  scale_size(range = c(0.25, 5)) +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  guides(colour = guide_legend(order = 1, nrow = 2), 
         size = guide_legend(order = 2, nrow = 3),
         shape = guide_legend(order = 3, nrow = 3)) +
  theme(legend.position = "bottom",
        legend.justification = c("left", "center")) +
  main_theme 

p1_2 <- ggplot(all, aes(POS)) + 
  geom_density(aes(color = FTYPE, y = ..scaled..)) +
  scale_color_manual(values = c("gold", "gray47")) +
  scale_y_continuous(position = "right") +
  xlab("") + ylab("Density of variants") +
  xlim(c(0, 5e+06)) +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  theme(legend.position = "none",
        axis.text.x=element_blank()) +
  main_theme

aligned_plots <- align_plots(p1_2, p1_1, align="hv", axis="tblr")
p_a <- ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])
###############################################################################
#### Panel B, variants between 1912 - 1945 Kb
### get region from peg.1535 to peg.1563
library(gggenes) 

gff <-read.delim("../00.data/annot_RAST_simp.gff", header = F, sep = "\t")
colnames(gff) <- c("Start", "End", "Strand", "Gene_ID", "Annotation")

region <- all[(all$POS > 1910551) & (all$POS < 1944637), ]
region <- merge(region, map)

gene_lst <- paste0("peg.", seq(1535, 1563))
region_gff <- gff[gff$Gene_ID %in% gene_lst, ]

## Get the coordinate of this gene
this_start <- min(region_gff$Start)
this_end <- max(region_gff$End)

colnames(region_gff) <- c("start", "end", "strand", "gene", "protein")

region_gff$orientation <- ifelse(region_gff$strand == "+", 1, 0)
region_gff$molecule <- ""

p_b1 <- ggplot(region_gff, aes(xmin = start, xmax = end,
                             y = molecule, fill = protein, 
                             forward = orientation)) +
  geom_gene_arrow() + labs(y = "") +
  xlim(this_start, this_end) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  theme(legend.position = "none")

##AF of each variants from each compounds
p_b2 <- ggplot(region, aes(POS, Ratio)) +
  geom_point(aes(color = Concentration, shape = Effect_type,
                 size = Ratio), alpha = 0.8) +
  xlim(this_start, this_end) +
  scale_shape_manual(limits = map$Effect_type,
                     values = map$Shape, guide = F) +
  scale_color_manual(values = c("0" = "gray",
                                "25" = "#a0cbe8", 
                                "250" = "navyblue",
                                "50" = "#1170aa", 
                                "500" = "navyblue")) +
  
  scale_size(range = c(0.25, 5), guide = F) +
  labs(x = "", y = "Allele Frequence (AF)") +
  main_theme +
  theme(legend.position = "top", 
        legend.justification = c("left", "center")) +
  #guides(shape = guide_legend(nrow = 1)) +
  theme(legend.background = element_blank())

p_b <- plot_grid(p_b2, p_b1, nrow = 2, rel_heights = c(4, 1),
          align = "v", axis = "l")

fig2_ab <- plot_grid(p_a, p_b, nrow = 2, align = "v", axis = "l", 
                  rel_heights = c(1, 1.2), 
                  labels = c('a', "b"))
###############################################################################
#### Panel C
pos_gene <- read.table("../00.data/vcf_20231025_position_gene_map.txt",
                       header = F, sep = "\t")
colnames(pos_gene) <- c("POS", "Type", "Gene")
all <- merge(all, pos_gene)

this_gene_ID <- c("peg.961", "peg.962", "peg.963")

this_gene <- all[all$Gene %in% this_gene_ID, ]
this_gene <- merge(this_gene, map)

this_gff <- gff[gff$Gene_ID %in% this_gene_ID, ]

## Get the coordinate of this gene
this_start <- min(gff[gff$Gene_ID %in% this_gene_ID, ]$Start)
this_end <- max(gff[gff$Gene_ID %in% this_gene_ID, ]$End)

p1 <- ggplot(this_gene, aes(POS, Ratio)) +
  geom_point(aes(color = Compound, shape = Effect_type, size = Concentration),
             alpha = 0.9) + 
  scale_size_manual(values = c("0" = 1,
                               "25" = 1.8, 
                               "250" = 3.6,
                               "50" = 2.4, 
                               "500" = 3.6)) +
  scale_color_manual(values = c("PFNA" = "lightblue4",
                                "PFOA" = "yellow4",
                                "Others" = "gray")) +
  labs(x = "Position", y = "AF") +
  xlim(this_start, this_end) +
  guides(color = guide_legend(nrow = 2), size = guide_legend(nrow = 2)) +
  scale_shape_manual(limits = map$Effect_type, values = map$Shape, guide = F) +
  theme(legend.position = "top") +
  main_theme +
  theme(#legend.key = element_rect(fill = "transparent"),
        legend.background = element_blank())

colnames(this_gff) <- c("start", "end", "strand", "gene", "protein")
this_gff$orientation <- ifelse(this_gff$strand == "+", 1, 0)
this_gff$molecule <- " "
p2 <- ggplot(this_gff, aes(xmin = start, xmax = end,
                           y = molecule, fill = protein, 
                           forward = orientation)) +
  geom_gene_arrow() + labs(y = "") +
  scale_fill_brewer(palette = "Greens") +
  theme_genes() +
  theme(legend.position = "none")

p_c <- plot_grid(p1, p2, nrow = 2, align = "v",  axis = "l", 
                 rel_heights = c(4, 1))

###############################################################################
#### Panel D
this_gene_ID <- c( "peg.1416", "peg.1417")

this_gene <- all[all$Gene %in% this_gene_ID, ]
this_gene <- merge(this_gene, map)

## Get the coordinate of this gene
this_gff <- gff[gff$Gene_ID %in% this_gene$Gene, ]
this_start <- min(this_gff$End)
this_end <- max(this_gff$End)

p1 <- ggplot(this_gene, aes(POS, Ratio)) +
  geom_point(aes(color = Compound, shape = Effect_type, size = Concentration),
             alpha = 0.9) + 
  xlim(this_start, this_end) +
  labs(x = "", y = "AF") +
  scale_size_manual(values = c("0" = 1,
                               "25" = 1.8, 
                               "50" = 2.4), guide = F) +
  scale_color_manual(values = c("Xanthan_gum" = "slateblue3",
                                "Others" = "gray")) +
  scale_shape_manual(limits = map$Effect_type, values = map$Shape, guide = F) +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 2)) +
  main_theme +
  theme(legend.background = element_blank())

colnames(this_gff) <- c("start", "end", "strand", "gene", "protein")
this_gff$orientation <- ifelse(this_gff$strand == "+", 1, 0)
this_gff$molecule <- " "

p2 <- ggplot(this_gff, aes(xmin = start, xmax = end,
                           y = molecule, fill = protein, 
                           forward = orientation)) +
  xlim(this_start, this_end) +
  geom_gene_arrow() + labs(y = "") +
  scale_fill_brewer(palette = "Purples") +
  theme_genes() +
  theme(legend.position = "none",
        panel.background = element_blank(),
        plot.background = element_blank()) 

p_d <- plot_grid(p1, p2, nrow = 2, align = "v",  axis = "l", 
                 rel_heights = c(4, 1))

###############################################################################
#### Panel E
this_gene_ID <- c( "peg.1703", "peg.1704")

this_gene <- all[all$Gene %in% this_gene_ID, ]
this_gene <- merge(this_gene, map)

## Get the coordinate of this gene
this_gff <- gff[gff$Gene_ID %in% this_gene$Gene, ]
this_start <- min(this_gff$Start)
this_end <- max(this_gff$End)

p1 <- ggplot(this_gene, aes(POS, Ratio)) +
  geom_point(aes(color = Compound, shape = Effect_type, size = Concentration),
             alpha = 0.9) + 
  xlim(this_start, this_end) +
  labs(x = "", y = "AF") +
  scale_size_manual(values = c("0" = 1,
                               "25" = 1.8, 
                               "50" = 2.4), guide = F) +
  scale_color_manual(values = c("Loperamide" = "#ba7e45",
                                "Others" = "gray")) +
  scale_shape_manual(limits = map$Effect_type, values = map$Shape, guide = F) +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 2)) +
  main_theme +
  theme(legend.background = element_blank())

colnames(this_gff) <- c("start", "end", "strand", "gene", "protein")
this_gff$orientation <- ifelse(this_gff$strand == "+", 1, 0)
this_gff$molecule <- " "

p2 <- ggplot(this_gff, aes(xmin = start, xmax = end,
                           y = molecule, fill = protein, 
                           forward = orientation)) +
  xlim(this_start, this_end) +
  geom_gene_arrow() + labs(y = "") +
  scale_fill_brewer(palette = "Oranges") +
  theme_genes() +
  theme(legend.position = "none",
        panel.background = element_blank(),
        plot.background = element_blank()) 

p_e <- plot_grid(p1, p2, nrow = 2, align = "v",  axis = "l", 
                 rel_heights = c(4, 1))

###############################################################################
####
fig2_cde <- plot_grid(p_c, p_d, p_e, nrow = 1, align = "h", axis = "b", 
                   labels = c("c", "d", "e"))

fig2_ab <- plot_grid(p_a, p_b, nrow = 2, align = "v", axis = "l", 
                     rel_heights = c(1, 1.2), 
                     labels = c('a', "b"))

fig2 <- plot_grid(fig2_ab, fig2_cde, nrow = 2, rel_heights = c(1.6, 1))

ggsave("../04.figures/Figure_2.pdf", fig2, width = 10, height = 9)
