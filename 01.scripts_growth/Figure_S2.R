library(curatedMetagenomicData)  # v3.8.0
library(dplyr)
library(stringr)
source("plot_settings.R")

adult_stool_h <-
  filter(sampleMetadata, age >= 18) |>
  filter(body_site == "stool") |>
  filter(disease == "healthy") |>
  #select(where(~ !all(is.na(.x)))) |>
  returnSamples("relative_abundance")

## Relative abundance matrix
ra_mat <- assay(adult_stool_h)
dim(ra_mat)  # 1478 5538

ra_mat <- as.matrix(apply(ra_mat, 2, function(x) x/sum(x)))

ra_occu <- data.frame(Taxonomy = rownames(ra_mat),
                      RA_mean = rowMeans(ra_mat),
                      RA_median = rowMedians(ra_mat),
                      Prevalence = rowSums(ra_mat>0)/ncol(ra_mat))

tax_levels <- str_split_fixed(ra_occu$Taxonomy, "\\|", 7)
colnames(tax_levels) <- c("Kingdom", "Phylum", "Class", "Order", "Family",
                          "Genus", "Species")
tax_levels <- as.data.frame(tax_levels)
tax_levels$Taxonomy <- ra_occu$Taxonomy

ra_occu <- merge(ra_occu, tax_levels)
ra_occu <- ra_occu[order(ra_occu$RA_mean, decreasing = T), ]
ra_occu_highlight <- ra_occu[1:5, ]

p_s2 <- ggplot(ra_occu, aes(RA_mean, Prevalence, color = Phylum)) + 
  geom_point() + 
  geom_point(data = ra_occu_highlight, size = 3) +
  geom_text(data = ra_occu_highlight,
            aes(label = Species),
            fontface = "italic",
            check_overlap = F, angle = 0, 
            nudge_x = -0.001, nudge_y = 0.03) +
  main_theme + 
  labs(x = "Average relative abundance") + 
  theme(legend.position = "bottom") +
  scale_color_manual(values = getOI(length(unique(ra_occu$Phylum)))) +
  guides(col = guide_legend(nrow = 3))

ggsave("../04.figures/Figure_S2.pdf", p_s2, height = 6, width = 13)
