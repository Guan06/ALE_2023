source("settings.R")

###############################################################################
# read in all variants called by freebayes
all <- readRDS("../04.results/20250910_all_variants_redo_qual.rds")

### Assign Plate as 2 for Metformin, as its DNA was added in different plates for sequencing, 
### during ALE experiment it was in plate 2, same for DMSO controls in column 9 and well A7

all[all$Compound == "Metformin", ]$Plate <- 2
all[(all$Compound == "DMSO_control" & all$Column == 9), ]$Plate <- 2
all[(all$Compound == "DMSO_control" & all$Well == "A7"), ]$Plate <- 2


#### Filter the variants with supporting reads < 5
all <- all[all$Alt_depth >= 5, ]
#### Filter the variants detected in parental samples 
all <- all[all$Ratio >= 0.05, ]

#### Identify parental variants 
parental <- all[all$Compound == "NT5002", ]

stat_prevalence <- parental %>% group_by(POS) %>% 
  summarise(Prevalence = length(unique(Sample_ID)))  

parental <- merge(parental, stat_prevalence)

parental$Prevalence <- as.character(parental$Prevalence)

ggplot(parental, aes(Sample_ID, Ratio, size = log(QUAL+1))) +
  geom_point(aes(color = Sample_ID, shape = Prevalence), position = "jitter") + 
  theme_bw() + labs(y = "AF") + scale_shape_manual(values = c(1, 2))

### For the identification of parental genetic variants, we applied following rules
### after manually check the mapping with IGV:
# - all variants with AF > 0.5. 
# - variants with prevalence == 2 (in triangle in the plot above). 
# - variants with prevalence == 1 BUT good quality (Q > 10). 
# - few exceptions confirmed that should be an existing variant but does not fulfill
# above criteria but we "see" them: 1,939,963; 4,229,448; and 4,452,307   

parental_vars <- parental[(parental$Ratio >= 0.5) | 
                            (parental$Prevalence == "2") | 
                            ((parental$Prevalence == "1") & (parental$QUAL >= 10)) |
                            (parental$POS %in% c(1939963, 4229448, 4452307)), ]

parental_vars <- unique(parental_vars$POS)
parental_vars
# cat(parental_vars, sep = ", ")
### 16 positions were identified as parental variants

### Remove parental variants from evolved populations
evol_pop <- all[all$Compound != "NT5002", ]
evol_pop <- evol_pop[!(evol_pop$POS %in% parental_vars), ]
## 1173

### Filter out positions detected in > 3 compounds
evol_stat <- evol_pop %>% group_by(POS) %>% 
  summarise(Comp_Prevalence = n_distinct(Compound))

evol_stat_multi <- unique(evol_stat[evol_stat$Comp_Prevalence >= 3, ]$POS)
#evol_stat_multi

evol_pop <- evol_pop[!(evol_pop$POS %in% evol_stat_multi), ]

#### Numbers used in the paper
print("Number of mutation")
print(dim(evol_pop)[1])

print("Number of unique mutation positions")
print(length(unique(evol_pop$POS)))

print("Median number of mutation")
print(median((evol_pop %>% group_by(Sample_ID) %>% summarise(n = n_distinct(POS)))$n))

saveRDS(evol_pop, "../04.results/20250910_evolved_population_variants_final.rds")
write.table(evol_pop, "../04.results/20250910_evolved_population_variants_final.txt",
            sep = "\t", row.names = F, quote = F)

## Distribution of AF in all mutations
ggplot(evol_pop, aes(Ratio)) + 
  geom_histogram(binwidth = 0.01, fill = NA, color = "gray57") + theme_bw()
## 7 variants with AF > 0.5; 5 of them are Xanthan gum - AcrR related variants

ggplot(evol_pop, aes(log(QUAL))) + 
  geom_histogram(binwidth = 1, fill = NA, color = "gray57") + theme_bw() #+
  # geom_vline(xintercept = log(10), color = "salmon", linetype = "dashed") +
  # geom_vline(xintercept = log(20), color =  "salmon", linetype = "dashed") +
  # geom_vline(xintercept = log(30), color =  "salmon", linetype = "dashed")
###############################################################################
### Figure 2a 
t1 <- evol_pop %>% group_by(POS) %>% summarise(Prevalence =
                                            length(unique(Sample_ID)),
                                          Average_AF = mean(Ratio))
t1_des <- unique(evol_pop[, c("POS", "FTYPE", "Effect")])
t1 <- merge(t1, t1_des)

map <- read.table("../00.data/vcf_20251021_effect_map.txt", 
                  header = T, sep = "\t")
map$Effect_type <- as.factor(map$Effect_type)
t1 <- merge(t1, map)

p1_1 <- ggplot(t1, aes(POS, Prevalence)) +
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=1771711, xmax=1772529, 
           alpha=0.8, fill="slateblue3") + 
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=2104963, xmax=2107776, 
           alpha=0.8, fill="#98dd94") + 
  geom_point(aes(color = FTYPE, size = Average_AF, shape = Effect_type),
             alpha = 0.8) + 
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

p1_2 <- ggplot(evol_pop, aes(POS)) + 
  geom_density(aes(color = FTYPE, 
                   y = (..scaled..))) + 
                   #y = after_stat(scaled))) +
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

p_a

library(gggenes) 

gff <-read.delim("../00.data/annot_RAST_simp.gff", header = F, sep = "\t")
colnames(gff) <- c("Start", "End", "Strand", "Gene_ID", "Annotation")

pos_gene <- read.table("../00.data/vcf_20231025_position_gene_map.txt",
                       header = F, sep = "\t")
colnames(pos_gene) <- c("POS", "Type", "Gene")
evol_pop <- merge(evol_pop, pos_gene, by = "POS", all.x = T)

#### Panel B - AcrR 
this_gene_ID <- c("peg.1416", "peg.1417")

this_gene <- evol_pop[evol_pop$Gene %in% this_gene_ID, ]
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
  scale_shape_manual(limits = map$Effect_type, values = map$Shape, guide = "none") +
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

p_b <- plot_grid(p1, p2, nrow = 2, align = "v",  axis = "l", 
                 rel_heights = c(4, 1))

###############################################################################
#### Panel C
this_gene_ID <- c( "peg.1703", "peg.1704")

this_gene <- evol_pop[evol_pop$Gene %in% this_gene_ID, ]
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
  scale_color_manual(values = c("Loperamide" = "#98dd94",
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
  scale_fill_manual(values = c("#98dd94", "#72a66f")) +
  theme_genes() +
  theme(legend.position = "none",
        panel.background = element_blank(),
        plot.background = element_blank()) 

p_c <- plot_grid(p1, p2, nrow = 2, align = "v",  axis = "l", 
                 rel_heights = c(4, 1))
p_c
###############################################################################
####
p_bc <- plot_grid(p_b, p_c, nrow = 1, align = "h", axis = "b", 
                  labels = c("b", "c"))

fig2 <- plot_grid(p_a, p_bc, nrow = 2, align="v", axis="lr", labels = c("a", ""),
                  rel_heights = c(1, 1.2))

fig2

ggsave("../05.figures/Figure_2.pdf", fig2, width = 7, height = 6)
