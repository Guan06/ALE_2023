## This scritps identify the "parental" variants in ALE 1.0 and 2.0
## and decide if they should be kept
###############################################################################
all <- readRDS("../00.data/vcf_20230829_all_variants_AF_Ratio_ALE1.rds")
all$Plate <- as.character(all$Plate)
parental <- all[all$Compound == "NT5002", ]

### Plot the position and prevalence of variants with AF > 0.05
af_02 <- parental[parental$Ratio > 0.05, ]
t1 <- af_02 %>% group_by(POS) %>% summarise(Prevalence =
                                              length(unique(Sample_ID)),
                                            Average_AF = mean(Ratio))

t1_des <- unique(af_02[, c("POS", "FTYPE")])
t1 <- merge(t1, t1_des)
colnames(t1)[1] <- "Position"

p1_1 <- ggplot(t1, aes(Position, Prevalence)) +
  geom_point(aes(color = FTYPE, size = Average_AF),
             alpha = 0.7, shape = 16) + 
  scale_color_manual(values = c("gold", "gray47")) +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  theme(legend.position = "top", legend.justification = c("left", "center")) +
  #  xlim(c(0, 5e+06)) +
  main_theme 

p1_2 <- ggplot(af_02, aes(POS)) + 
  geom_density(aes(color = FTYPE)) +
  scale_color_manual(values = c("gold", "gray47")) +
  scale_y_continuous(position = "right") +
  xlab("") + ylab("Density of prevalence") +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  theme(legend.position = "none") +
  main_theme

aligned_plots <- align_plots(p1_2, p1_1, align="hv", axis="tblr")
p1 <- ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])


# p1 <- ggplot(t1, aes(Position, Prevalence)) +
#   geom_point(aes(color = FTYPE, size = Average_AF),
#              alpha = 0.7, shape = 16) + 
#   scale_color_manual(values = c("gold", "gray47")) +
#   #scale_size_manual(values = c(1, 2)) +
#   theme(legend.position = "top") +
#   xlim(c(0, 5e+06)) +
#   main_theme 

ggsave("../04.figures/Appendix_fig_01_parental_variants_AF.pdf", p1, 
       width = 7, height = 3)
###############################################################################
parental_lst1 <- t1[(t1$Prevalence == 2) & (t1$Average_AF > 0.9), ]$Position

parental_sub <- all[all$POS %in% parental_lst1, ]
parental_sub <- parental_sub[!is.na(parental_sub$Solvent), ]
p2 <- ggplot(parental_sub, aes(interaction(Plate, Compound), fill = Concentration)) + 
    geom_bar(stat = "count", position = position_dodge2(preserve = "single")) + 
    main_theme +  
    scale_fill_manual(values = c("0" = "gray",
                                "25" = "#a0cbe8", 
                                "50" = "#1170aa")) +
    facet_grid(rows = vars(POS), cols = vars(Solvent),
               scales = "free", space = "free") +
    geom_hline(yintercept = 4, color = "gold", linetype = "dashed") +
    geom_hline(yintercept = 2, color = "salmon", linetype = "dashed") +
    labs(y = "Number of samples", x = "") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

ggsave("../04.figures/Appendix_fig_02_parental_variants_in_ALE1.pdf", p2, 
       width = 11, height = 8)

###############################################################################
ale2 <- readRDS("../00.data/vcf_20230814_all_variants_ALE2.rds")

af_005 <- ale2[ale2$Ratio > 0.05, ]
t2 <- af_005 %>% group_by(POS)  %>% summarise(Prevalence =
                                        length(unique(Sample_ID)),
                                      Average_AF = mean(Ratio))
t2_des <- unique(af_005[, c("POS", "FTYPE")])
t2 <- merge(t2, t2_des)
colnames(t2)[1] <- "Position"

p3_1 <- ggplot(t2, aes(Position, Prevalence)) +
  geom_point(aes(color = FTYPE, size = Average_AF),
             alpha = 0.7, shape = 16) + 
  scale_color_manual(values = c("gold", "gray47")) +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  theme(legend.position = "top") +
#  xlim(c(0, 5e+06)) +
  main_theme 

p3_2 <- ggplot(af_005, aes(POS)) + 
  geom_density(aes(color = FTYPE)) +
  scale_color_manual(values = c("gold", "gray47")) +
  scale_y_continuous(position = "right") +
  xlab("") + ylab("Density of prevalence") +
  theme_minimal_hgrid(color = "gray88", font_size = 10) +
  theme(legend.position = "none") +
  main_theme

aligned_plots <- align_plots(p3_2, p3_1, align="hv", axis="tblr")
p3 <- ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])

ggsave("../04.figures/Appendix_fig_03_common_variants_AF_ALE2.pdf", p3, 
       width = 7, height = 3)

###############################################################################
com_lst <- t2[(t2$Prevalence == 12) & (t2$Average_AF > 0.9), ]$Position
ale2_com <- ale2[ale2$POS %in% com_lst, ]
write.table(ale2_com,
            "../05.tables/Appendix_tab_01_common_variants_ALE2.txt",
            quote = F, sep = "\t", row.names = F)
