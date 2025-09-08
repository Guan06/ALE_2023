source("../02.scripts_variants/settings.R")
library(ggridges)

od_iso <- read.table("../00.data//20240912_iso_OD_parental_and_XG.txt",
                     header = T, sep = "\t")
od_iso$Date <- as.character(od_iso$Date)

### Remove obvious outlier caused by experimental error? e.g. the last row 
### (DMSO) was not inoculated, therefore no growth there
od_iso<- od_iso[!(od_iso$Row == "H" & od_iso$OD < 0.4), ]
od_iso <- od_iso[!(od_iso$Row == "D" & od_iso$OD > 0.4 & od_iso$Antibiotic == "Metronidazole"), ]

t1 <- od_iso %>% group_by(Group, Antibiotic) %>% 
  mutate(Min_OD = min(OD), Max_OD = max(OD)) %>% 
  mutate(MIC80_OD = (Max_OD - Min_OD) * 0.2 + Min_OD)


t1 <- t1[t1$OD > t1$MIC80_OD, ]

t1 <- t1 %>% group_by(Group, Antibiotic) %>% mutate(MIC80 = max(Concentration))

t2 <- unique(t1[, c("Group", "Compound", "Compound_conc2", 
                    "MIC80", "Antibiotic", "Date", "Population")])

t2 <- t2 %>% mutate(Rep = case_when(
  Date == "20240621" ~ 'Parental',
  (Date == "20240720" | Date == "20240705") ~ 'Batch 1',
  (Date == "20240802" | Date == "20240808") ~ 'Batch 2',
))

p_3b <- ggplot(t2, aes(y = interaction(Compound_conc2, Rep), x = MIC80,
               fill = Compound_conc2)) +
  geom_jitter(aes(color = Compound_conc2, shape = Population), alpha = 0.6,
              width = 0, height = 0.3) +
  geom_density_ridges(scale = 0.5, quantile_lines = F, quantiles = 2,
                      color = "gray67", lwd = 0.3, alpha = 0.7) +
  facet_grid(cols=vars(Antibiotic), rows=vars(Rep),
             scales = "free", space = "free_y") +
  scale_color_manual("Population (compound and concentration)", values = cc_color) +
  scale_fill_manual(values = cc_color, guide = F) +
  scale_shape_manual(values = cc_shape) +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape = guide_legend(order = 2, nrow = 2)) +
  #  scale_shape_manual(values = c(5, 2, 1)) +
  main_theme +
  theme(axis.text.y = element_blank(),
        legend.position = "top") +
  labs(y = "Compound and concentration", x = "MIC")

###############################################################################
## 2025.06.10 manual check the outlier data points
od_dox <- od_iso[od_iso$Group == "iso418-3_C2_20240808_8", ]
cutoff <- (max(od_dox$OD) - min(od_dox$OD)) * 0.2 + min(od_iso$OD)

ggplot(od_dox, aes(Row, OD)) +
  geom_point() + geom_line() + 
  geom_hline(yintercept = cutoff, color = "salmon") +
  main_theme

## 2025.06.30 manual check the outlier data points
## xanthan gum 25 with MIC of 200
od_dox <- od_iso[od_iso$Group == "iso418-3_E12_20240808_2", ]
cutoff <- (max(od_dox$OD) - min(od_dox$OD)) * 0.2 + min(od_iso$OD)

ggplot(od_dox, aes(Row, OD)) +
  geom_point() + geom_line() + 
  geom_hline(yintercept = cutoff, color = "salmon") +
  main_theme

### The cutoff is underestimated due to the OD read in row E
### MIC should be the concentration of row E (0.125)
t2[t2$Group == "iso418-3_C2_20240808_8", ]$MIC80 <- 0.125

#####

t2_total <- t2 %>% group_by(Compound_conc2, Rep, Antibiotic) %>%
  summarise(Num = n())

t2_mic <- t2 %>% 
  group_by(Compound_conc2, Rep, MIC80, Antibiotic) %>%
  summarise(MIC_num = n())

t2_mic <- merge(t2_mic, t2_total, 
                by = c("Antibiotic", "Compound_conc2", "Rep"))

t2_mic$Ratio <- t2_mic$MIC_num / t2_mic$Num 
t2_mic$logX <- log2(t2_mic$MIC80)

t2_mic_median <- t2 %>% 
  filter(Rep == "Parental") %>% 
  group_by(Antibiotic) %>%
  summarise(Median = median(MIC80))

t2_mic <- t2_mic %>% mutate(Rep = case_when(
  Rep == "Parental" ~ 'Parental',
  (Rep == "Replicate 1") ~ 'Batch 1',
  (Rep == "Replicate 2") ~ 'Batch 2',
))

p_3b_v2 <- ggplot(t2_mic, aes(x = MIC80, y = interaction(Compound_conc2, Rep))) + 
  geom_point(aes(size = Ratio, color = Compound_conc2)) +  
  scale_size_continuous(range = c(0.5, 5)) +
  facet_grid(cols=vars(Antibiotic),  
             rows = vars(Rep), scales = "free", space = "free_y") + 
  geom_vline(data = t2_mic_median, aes(xintercept = Median), 
             linetype = "solid", color = "salmon", size = 1.8, alpha = 0.4) +
  scale_fill_manual(values = cc_color, guide = F) + 
  scale_color_manual("Population (compound and concentration)", 
                     values = cc_color) + 
  theme_bw() +
  #main_theme +
  theme(axis.text.y = element_blank(),
        legend.position = "top",
        panel.spacing.x = unit(0.8, "lines")) +
 
  labs(y = "", 
       x = "MIC (Minimium Inhibitory Concentration µg / mL)") +
  scale_x_continuous(expand = c(0.1, 0)) 

p_3b_v2
ggsave("../05.figures/Figure_3b_v2.pdf", width = 10, height = 5)

################################################################################
############################Plot each population 

t2_total <- t2 %>% group_by(Compound_conc2, Population, Rep, Antibiotic) %>%
  summarize(Num = n())

t2_mic <- t2 %>% 
  group_by(Compound_conc2, Population, Rep, MIC80, Antibiotic) %>%
  summarize(MIC_num = n())

t2_mic <- merge(t2_mic, t2_total, 
                by = c("Antibiotic", "Compound_conc2", "Population","Rep"))

t2_mic$Ratio <- t2_mic$MIC_num / t2_mic$Num 
t2_mic$logX <- log2(t2_mic$MIC80)

ggplot(t2_mic, aes(x = logX, y = Population)) + 
  geom_point(aes(size = Ratio, color = Compound_conc2)) +  
  scale_size_continuous(range = c(0.4, 6)) +
  #scale_x_continuous(breaks = log2(pretty(2 ** (t2_mic$logX))), 
  #                   labels = pretty(2 ** (t2_mic$logX))) +
  facet_grid(cols=vars(Antibiotic),  
             rows = vars(Compound_conc2), scales = "free", space = "free_y") + 
  #facet_wrap(vars(Rep, Antibiotic), ncol = 6, scales = "free",
  #           strip.position = "top") +
  #scale_x_continuous(breaks = log2(pretty(t2_mic$MIC80)),
  #                   labels = pretty(t2_mic$MIC80)) +
  scale_fill_manual(values = cc_color, guide = F) + 
  scale_color_manual("Population (compound and concentration)", 
                     values = cc_color) + 
  theme_bw() +
  theme(axis.text.y = element_blank(),
        legend.position = "top",
        panel.spacing.x = unit(0.8, "lines")) +
  labs(y = "Compound and concentration", x = "log2(MIC)") #+
#scale_x_continuous(expand = c(-0.2, 0)) +

#ggsave("../05.figures/Figure_3b_v3.pdf", width = 10, height = 6)
################################################################################

### Supplementary Tables S4 - S6
###############################################################################
## Amoxicillin
amo <- t2[t2$Antibiotic == "Amoxicillin", ]
amo_stat <- amo %>% group_by(Compound, Population, MIC80) %>% summarise(n = n())
amo_stat2 <- amo_stat %>% pivot_wider(names_from = MIC80, values_from = n)
print(amo_stat2)

ery <- t2[t2$Antibiotic == "Erythromycin", ]
ery_stat <- ery %>% group_by(Compound, Population, MIC80) %>% summarise(n = n())
ery_stat2 <- ery_stat %>% pivot_wider(names_from = MIC80, values_from = n)
ery_stat2[is.na(ery_stat2)] <- 0
print(ery_stat2)

## Metronidazole
met <- t2[t2$Antibiotic == "Metronidazole", ]
met_stat <- met %>% group_by(Compound, Population, MIC80) %>% summarise(n = n())
met_stat2 <- met_stat %>% pivot_wider(names_from = MIC80, values_from = n)
met_stat2[is.na(met_stat2)] <- 0
print(met_stat2)



### Previous Figure 3c
###############################################################################
library(pracma)
auc<- od_iso %>% group_by(Group, Antibiotic, Population) %>% 
  arrange(Concentration) %>% 
  summarise(AUC = trapz(Concentration, OD))

des <- unique(od_iso[, colnames(od_iso) %in% c("Group", "Date",
                                               "Antibiotic", "Compound_conc2")])

auc_des <- merge(auc, des)

shapes <- c("Parental" = 1, 
            "Plate4_A4" = 1, "Plate4_B4" = 2, "Plate4_C4" = 5, "Plate4_D4" = 8,
            "Plate4_E4" = 1, "Plate4_F4" = 2, "Plate4_G4" = 5, "Plate4_H4" = 8,
            "Plate4_A5" = 1, "Plate4_B5" = 2, "Plate4_C5" = 5, "Plate4_D5" = 8)

p_3c <- ggplot(auc_des, aes(Compound_conc2, AUC)) +
  geom_jitter(aes(shape = Population, color = Compound_conc2), alpha = 0.5) +
  geom_boxplot(aes(color = Compound_conc2), fill = NA, 
               alpha = 0.8, outlier.shape = NA) +
  scale_shape_manual(values = shapes) +
  facet_wrap(~Antibiotic, scales = "free_y", nrow = 1) +
  #theme_bw() +
  main_theme +
  labs(x = "", y = "IRE") +
  scale_color_manual(values = cc_color) +
  theme(axis.text.x = element_blank(),
        legend.position = "none")
p_3c
#ggsave("../05.figures/Figure_3c_v2.pdf", width = 10, height = 3)

p_3bc <- plot_grid(p_3b_v2, p_3c, ncol = 1, rel_heights = c(2.4, 1), 
                   labels = c('b', 'c'), align = "v", axis = "lr")

p_3bc
ggsave("../05.figures/Figure_3bc.pdf", width = 10, height = 8)

###############################################################################
### Significance test, results summarised in Table S7
abx_lst <- unique(auc_des$Antibiotic)

for (i in abx_lst) {
  this <- auc_des[auc_des$Antibiotic == i, ]
  print(i)
  print(pairwise.wilcox.test(this$AUC, this$Compound_conc2, 
                             p.adjust.method = "bonf"))
}

###############################################################################
### Figure 3c - v2
###############################################################################
library(pracma)
auc<- od_iso %>% group_by(Group, Antibiotic, Population) %>% 
  arrange(Concentration) %>% 
  summarise(AUC = trapz(Concentration, OD))

des <- unique(od_iso[, colnames(od_iso) %in% c("Group", "Date",
                                               "Antibiotic", "Compound_conc2")])

auc_des <- merge(auc, des)
auc_des <- auc_des %>% mutate(Rep = case_when(
  Date == "20240621" ~ 'Parental',
  (Date == "20240720" | Date == "20240705") ~ 'Batch 1',
  (Date == "20240802" | Date == "20240808") ~ 'Batch 2',
))

p_3c_v2 <- ggplot(auc_des, aes(Population, AUC)) +
  geom_jitter(aes(shape = Rep, color = Compound_conc2), alpha = 0.5) +
  geom_boxplot(aes(color = Compound_conc2), fill = NA, 
               alpha = 0.8, outlier.shape = NA) +
  #scale_shape_manual(values = shapes) +
  facet_wrap(~Antibiotic, scales = "free_y", nrow = 2) +
  scale_x_discrete(limits=c("Parental",
                            "Plate4_A5", "Plate4_B5", "Plate4_C5",
                            "Plate4_A4", "Plate4_B4", "Plate4_C4", "Plate4_D4",
                            "Plate4_E4", "Plate4_F4", "Plate4_G4", "Plate4_H4")) +
  main_theme +
  scale_shape_manual(values = c(2, 1, 5)) +
  labs(x = "", y = "IRE") +
  scale_color_manual(values = cc_color) +
  theme(axis.text.x = element_blank(),
        legend.position = "top")


p_3bc_v2 <- plot_grid(p_3b_v2, p_3c_v2, ncol = 1, rel_heights = c(1.15, 1), 
                   labels = c('b', 'c'), align = "v", axis = "lr")

p_3bc_v2
ggsave("../05.figures/Figure_3bc_v2.pdf", width = 10, height = 10)

abx_lst <- unique(auc_des$Antibiotic)

cmp <- c()
for (i in abx_lst) {
  this <- auc_des[auc_des$Antibiotic == i, ]
  print(i)
  print(pairwise.wilcox.test(this$AUC, this$Population, 
                             p.adjust.method = "bonf"))
  this_sig <- pairwise.wilcox.test(this$AUC, this$Population, 
                                   p.adjust.method = "bonf")$p.value
  this_sig <- as.matrix(this_sig)
  
  this_tab <- data.frame(as.table(this_sig))[lower.tri(this_sig, diag = TRUE), ]
  colnames(this_tab) <- c("Pop_1", "Pop_2", "P_adj")
  this_tab$Antibiotic <- i
  
  cmp <- rbind(cmp, this_tab)
}

cmp$Sig <- ifelse(cmp$P_adj < 0.05, "Yes", "No")
write.table(cmp, "../04.results/Table_EV10_full.txt", quote = F, row.names = F)

## filter out the comparison doesn't contain control group
ctr_lst <- c("Parental", "Plate4_A5", "Plate4_B5", "Plate4_C5")
cmp2 <- cmp[(cmp$Pop_1 %in% ctr_lst) |(cmp$Pop_2 %in% ctr_lst), ]
write.table(cmp2, "../04.results/Table_EV10.txt", quote = F, row.names = F)


################################################################################
# Function to assign groups based on significant differences 
source("function_sig_grp.R")
abx_lst <- unique(cmp$Antibiotic)
cmp_grp <- c()
for (abx in abx_lst){
  this_cmp <- cmp[cmp$Antibiotic == abx, ]
  
  #this_grp <- cld(this_cmp)
  this_grp <- assignLetters(this_cmp)
  
  
  this_grp$Antibiotic <- abx
  
  cmp_grp <- rbind(cmp_grp, this_grp)
}
write.table(cmp_grp, "../04.results/Table_EV10_group.txt",
            quote = F, row.names = F)
###############################################################################




###############################################################################
t <- auc_des[auc_des$Antibiotic == "Erythromycin" & auc_des$Population == "Plate4_H4", ]
## iso418-5_G3 behave strangely

od_dox <- od_iso[grepl("iso418-5_G3", od_iso$Group), ]
cutoff <- (max(od_dox$OD) - min(od_dox$OD)) * 0.2 + min(od_iso$OD)

t1 <- ggplot(od_dox, aes(Row, OD)) +
  geom_point() + geom_line(aes(group = Group)) + 
  facet_wrap(~Antibiotic, scales = "free_y", nrow = 1) +
  geom_hline(yintercept = cutoff, color = "salmon") +
  main_theme

od_dox <- od_iso[grepl("iso418-5_G1", od_iso$Group), ]
cutoff <- (max(od_dox$OD) - min(od_dox$OD)) * 0.2 + min(od_iso$OD)

t2 <- ggplot(od_dox, aes(Row, OD)) +
  geom_point() + geom_line(aes(group = Group)) + 
  facet_wrap(~Antibiotic, scales = "free_y", nrow = 1) +
  geom_hline(yintercept = cutoff, color = "salmon") +
  main_theme

plot_grid(t1, t2, nrow = 2)

tt <- auc_des[auc_des$Antibiotic == "Doxycycline" & auc_des$Population == "Plate4_H4", ]
## iso418-5_G3 behave strangely

od_dox <- od_iso[grepl("iso418-5_G3", od_iso$Group), ]
cutoff <- (max(od_dox$OD) - min(od_dox$OD)) * 0.2 + min(od_iso$OD)
###############################################################################