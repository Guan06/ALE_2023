## This is a potential supplementary of Figure 2C, where all the variants that 
## appeared in parental strains that have AF > 0.1 will be removed from this plot

#### Panel C, variants between 1912 - 1945 Kb
### get region from peg.1535 to peg.1563
library(gggenes) 

gff <-read.delim("../00.data/annot_RAST_simp.gff", header = F, sep = "\t")
colnames(gff) <- c("Start", "End", "Strand", "Gene_ID", "Annotation")

# read in all variants detected in ALE 1.0 and ALE 2.0 with AF > 0.05
all <- readRDS("../04.results/20231023_all_005.rds")
map <- read.table("../00.data/vcf_20231024_effect_map.txt", 
                  header = T, sep = "\t")
map$Effect_type <- as.factor(map$Effect_type)

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

## Get the list of position where variants were detected in parental strains

region_parental <- region[region$Sample_ID %in% c("Plate1B9", "Plate1C9"), ]

region_parental <- region_parental[region_parental$Ratio > 0.1, ]
lst <- unique(region_parental$POS)

region_filter <- region[!(region$POS %in% lst), ]

##AF of each variants from each compounds
p_c2 <- ggplot(region_filter, aes(POS, Ratio)) +
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

ggsave("../05.figures/Figure_2c_filter.pdf", p_c, width = 10, height = 3)

