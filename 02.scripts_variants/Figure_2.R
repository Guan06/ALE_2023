source("settings.R")

###############################################################################
# read in all variants detected in ALE 1.0 with AF > 0.05
all <- readRDS("../04.results/20231023_all_005.rds")
all <- all[all$Exp == "ALE1", ]
t1 <- all %>% group_by(POS) %>% summarise(Prevalence =
                                            length(unique(Sample_ID)),
                                          Average_AF = mean(Ratio))
t1_des <- unique(all[, c("POS", "FTYPE", "Effect")])
t1 <- merge(t1, t1_des)

map <- read.table("../00.data/vcf_20231024_effect_map.txt", 
                  header = T, sep = "\t")
map$Effect_type <- as.factor(map$Effect_type)
#map$Shape <- as.factor(map$Shape)

t1 <- merge(t1, map)

p1_1 <- ggplot(t1, aes(POS, Prevalence)) +
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=1910551, xmax=1944637, 
           alpha=0.3, fill="skyblue") + 
  geom_point(aes(color = FTYPE, size = Average_AF, shape = Effect_type),
             alpha = 0.7) + 
  scale_shape_manual(limits = map$Effect_type,
                     values = map$Shape) +
  scale_color_manual(values = c("gold", "gray47")) +
  xlim(c(0, 4.685e+06)) + xlab("") +
  scale_size(range = c(0.25, 5), limits = c(0.05, 1)) +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  guides(colour = guide_legend(order = 2, nrow = 1), 
         size = guide_legend(order = 3, nrow = 1),
         shape = guide_legend(order = 1, nrow = 1)) +
  main_theme + theme(legend.position = "top")

legend <- get_plot_component(p1_1, 'guide-box-top', return_all = TRUE)
p1_1 <- p1_1 + theme(legend.position = "none")

p1_2 <- ggplot(all, aes(POS)) + 
  geom_density(aes(color = FTYPE, y = ..scaled..)) +
  scale_color_manual(values = c("gold", "gray47")) +
  scale_y_continuous(position = "right") +
  xlab("") + ylab("Density of variants") +
  xlim(c(0, 4.685e+06)) +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  theme(legend.position = "none",
        axis.text.x=element_blank()) +
  main_theme

aligned_plots <- align_plots(p1_2, p1_1, align="hv", axis="tblr")
p_a <- ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])
p_a <- plot_grid(legend, p_a, nrow = 2, rel_heights = c(0.2, 1))

###############################################################################
# read in all variants detected in ALE 1.0 with AF > 0.5
all_05 <- readRDS("../04.results/20231023_all_05.rds")
all_05 <- all_05[all_05$Exp == "ALE1", ]
t1 <- all_05 %>% group_by(POS) %>% summarise(Prevalence =
                                            length(unique(Sample_ID)),
                                          Average_AF = mean(Ratio))
t1_des <- unique(all_05[, c("POS", "FTYPE", "Effect")])
t1 <- merge(t1, t1_des)
t1 <- merge(t1, map)

p2_1 <- ggplot(t1, aes(POS, Prevalence)) +
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=1910551, xmax=1944637, 
           alpha=0.3, fill="skyblue") + 
  geom_point(aes(color = FTYPE, size = Average_AF, shape = Effect_type),
             alpha = 0.7) + 
  scale_color_manual(values = c("gold", "gray47")) +
  xlab("Position") +
  scale_shape_manual(limits = map$Effect_type,
                     values = map$Shape) + 
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  scale_size(range = c(0.25, 5), limits = c(0.05, 1)) +
  xlim(c(0, 4.685e+06)) +
  guides(colour = guide_legend(order = 1, nrow = 1), 
         size = guide_legend(order = 2, nrow = 1),
         shape = guide_legend(order = 3, nrow = 1)) +
  theme(legend.position = "bottom",
        legend.justification = c("left", "center")) +
  main_theme 

legend <- get_plot_component(p2_1, 'guide-box-bottom', return_all = TRUE)
p2_1 <- p2_1 + theme(legend.position = "none")

p2_2 <- ggplot(all_05, aes(POS)) + 
  geom_density(aes(color = FTYPE, y = ..scaled..)) +
  scale_color_manual(values = c("gold", "gray47")) +
  scale_y_continuous(position = "right") +
  xlab("") + ylab("Density of variants") +
  xlim(c(0, 4.685e+06)) +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  theme(legend.position = "none",
        axis.text.x=element_blank()) +
  main_theme

aligned_plots <- align_plots(p2_2, p2_1, align="hv", axis="tblr")
p2 <- ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])
p_b <- plot_grid(legend, p2, nrow = 2, rel_heights = c(0.2, 1))

###############################################################################
#### Panel C, variants between 1912 - 1945 Kb
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

p_c1 <- ggplot(region_gff, aes(xmin = start, xmax = end,
                             y = molecule, fill = protein, 
                             forward = orientation)) +
  geom_gene_arrow() + labs(y = "") +
  xlim(this_start, this_end) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  theme(legend.position = "none")

region <- region[!(region$Sample_ID %in% c("Plate1B9",
                                           "Plate1C9")), ]

##AF of each variants from each compounds
p_c2 <- ggplot(region, aes(POS, Ratio)) +
  geom_point(aes(color = Concentration, shape = Effect_type), 
             size = 3, alpha = 0.8) +
  xlim(this_start, this_end) +
  scale_shape_manual(limits = map$Effect_type,
                     values = map$Shape, guide = F) +
  scale_color_manual(values = c("0" = "gray",
                                "25" = "#a0cbe8", 
                                "250" = "navyblue",
                                "50" = "#1170aa", 
                                "500" = "navyblue")) +
  #scale_size(range = c(0.25, 5), guide = F) +
  labs(x = "", y = "Allele Frequence (AF)") +
  theme_bw() +
  theme(legend.position = "top", 
        legend.justification = c("left", "center")) +
  #guides(shape = guide_legend(nrow = 1)) +
  theme(legend.background = element_blank()) 

p_c <- plot_grid(p_c2, p_c1, nrow = 2, rel_heights = c(4, 1),
          align = "v", axis = "l")

fig2_abc <- plot_grid(p_a, p_b, p_c, nrow = 3, align = "v", axis = "lr", 
                  rel_heights = c(1, 1, 1.3), 
                  labels = c('a', "b", "c"))

###############################################################################
###############################################################################
pos_gene <- read.table("../00.data/vcf_20231025_position_gene_map.txt",
                       header = F, sep = "\t")
colnames(pos_gene) <- c("POS", "Type", "Gene")
all <- merge(all, pos_gene)

#### Panel D - AcrR 
this_gene_ID <- c("peg.1416", "peg.1417")

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
fig2_de <- plot_grid(p_d, p_e, nrow = 1, align = "h", axis = "b", 
                   labels = c("d", "e"))


fig2 <- plot_grid(fig2_abc, fig2_de, nrow = 2, rel_heights = c(2.4, 1))

ggsave("../05.figures/Figure_2.pdf", fig2, width = 10, height = 12)


###############################################################################
#### Numbers used in the paper
all_noNT5002 <- all[all$Compound != "NT5002", ]
dim(all_noNT5002)
length(unique(all_noNT5002$POS))

number_of_var <- all_noNT5002 %>% group_by(Sample_ID) %>% summarise(n = n())
median(number_of_var$n)

number_of_samples <- all_noNT5002 %>% 
  group_by(Compound) %>% summarise(n = length(unique(Sample_ID)))

all_05_noNT5002 <- all_05[all_05$Compound != "NT5002", ]
dim(all_05_noNT5002)
length(unique(all_05_noNT5002$POS))

number_of_var_05 <- all_05_noNT5002 %>% group_by(Sample_ID) %>% summarise(n = n())
median(number_of_var_05$n)
