source("settings.R")

## Make sure to run pre_Figure_S8.R to get the  tri-nucleotides and
## their reverse compliment for all the position with detected variation.

all_snp <- read.delim2("../03.results/Figure_S8_all_snp.txt", header = T)
id_to_mut <- read.table("../00.data/20221127_ID_to_mutation_type.txt", 
                        header = T, sep = "\t")

all_snp_stat <- all_snp %>% group_by(Mutation, Sample_ID, .drop = F) %>% 
  summarise(Count = n())

###############################################################################
stat_mat <- all_snp_stat %>% pivot_wider(names_from = Sample_ID, values_from = Count)
stat_mat[is.na(stat_mat)] <- 0
###############################################################################

stat_AF <- all_snp%>% group_by(Compound, Mutation, .drop = T) %>%
  summarise(Count = n())


stat_AF <- merge(stat_AF, id_to_mut)
stat_AF$Substitution <- str_split_fixed(stat_AF$ID, "\\.", 2)[, 1]

comp_lst <- unique(stat_AF$Compound)

for (comp in comp_lst) {
  if (comp %in% c("Water_control", "DMSO_control")) next 
  this_comp <- stat_AF[stat_AF$Compound == comp, ]
  
  ## Get the solvent 
  solvent <- ifelse(sum(this_comp %in% c("Glyphosate", "Ethanol", "Xanthan_gum")),
                         "Water_control", "DMSO_control")
  this_solvent <- stat_AF[stat_AF$Compound == solvent, ]

  this_comp$Count_ratio <- this_comp$Count / sum(this_comp$Count)
  this_solvent$Count_ratio <- -(this_solvent$Count / sum(this_solvent$Count))
  
  this <- rbind(this_comp, this_solvent)
  this_p <- ggplot(this, aes(Mutation, Count_ratio)) +
    geom_bar(aes(fill = Substitution), 
             alpha = 0.7, color = "gray47",
             stat = 'identity', position = 'identity') +
    geom_hline(yintercept = 0) +
    labs(y = "Count %", x = "") +
    ggtitle(paste0(sum(this_comp$Count), 
                  " SNPs were detected in ", comp, " populations and ",
                  sum(this_solvent$Count),
                  " SNPs in ", solvent,  " populations.")) + 
    scale_fill_manual(values = palette_OkabeIto[1:6]) +
    scale_x_discrete(limits = id_to_mut$Mutation) +
    main_theme +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 5.4, 
                                     angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0("../04.figures/Figure_S8/", comp, ".pdf"),
           this_p, width = 9, height = 3)
}

###############################################################################
########## For PFNA and PFOA at different concentrations
pfas_snp <- all_snp[all_snp$Compound %in% c("PFNA", "PFOA", "DMSO_control"), ]

## Only get DMSO control from Plate3
t <- unique(pfas_snp$Sample_ID)[grepl("Plate", unique(pfas_snp$Sample_ID))]
t_other <- t[!grepl("Plate3", t)]
##

pfas_snp <- pfas_snp[!(pfas_snp$Sample_ID %in% t_other), ]

pfas_stat <- pfas_snp %>% group_by(Exp, Compound, Concentration, Mutation, ID,
                                   .drop = T) %>%
  summarise(Count = n())
  
pfas_stat <- pfas_stat %>% group_by(Compound, Concentration) %>% 
  mutate(Proportion = Count/sum(Count))
  

pfas_stat$Substitution <- str_split_fixed(pfas_stat$ID, "\\.", 2)[, 1]

p <- ggplot(pfas_stat, aes(Mutation, Proportion)) +
  geom_col(aes(fill = Substitution), alpha = 0.7, color = "NA") +
  #geom_point(color = "gray47", shape = 1) +
  scale_fill_manual(values = palette_OkabeIto[1:6]) +
  scale_x_discrete(limits = id_to_mut$Mutation) +
  #facet_grid(cols = vars(Compound), rows = vars(Concentration), scales = "free") +
  facet_wrap(~ Exp + Concentration + Compound, ncol = 2, scales = "free",
             labeller = label_wrap_gen(multi_line=FALSE)) +
  main_theme +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust=1))
ggsave("../04.figures/Figure_S8.pdf", p, width = 14, height = 8)
