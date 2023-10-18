## function sourced by Figure_S3.R, Figure_1.R
############################################################################
## significant test for DMSO vs Water control
get_sig_control <- function(x, group, compare) {
    lst <- as.character(unique(x[[group]]))
    len <- length(lst)
    sig <- c()
    for (i in 1 : (len - 1 )) {
        for (j in (i + 1) : len) {
            sig_i <- x[x[[group]] == lst[i], ]
            sig_j <- x[x[[group]] == lst[j], ]
            sig_ij <- wilcox.test(sig_i[[compare]], sig_j[[compare]],
                                  alternative = "two.sided")
            this_sig <- c(lst[i], lst[j], sig_ij$p.value)
            sig <- rbind(sig, this_sig)
        }
    }
    colnames(sig) <- c("Group1", "Group2", "Significance")
    rownames(sig) <- c()
    sig <- as.data.frame(sig)

    sig$Passage1 <- sapply((strsplit(sig$Group1, "_")), "[[", 1)
    sig$Passage2 <- sapply((strsplit(sig$Group2, "_")), "[[", 1)
    sig <- sig[sig$Passage1 == sig$Passage2, ]

    sig$FDR <- p.adjust(sig$Significance, method = "fdr")
    sig$FDR_sig <- ifelse(sig$FDR < 0.05, "Sig", "Non-sig")

    sig$Compare <- compare
    return(sig)
}

############################################################################
## significant test for compound vs control
get_sig <- function(x, group, compare) {
  lst <- as.character(unique(x[[group]]))
  len <- length(lst)
  sig <- c()
  for (i in 1 : (len - 1 )) {
    for (j in (i + 1) : len) {
      sig_i <- x[x[[group]] == lst[i], ]
      sig_j <- x[x[[group]] == lst[j], ]
      
      if (unique(sig_i$Compound) == unique(sig_j$Compound)) {next}
    
      ## check which one is control group
      if (unlist(strsplit(lst[j], "_"))[3] == "control"){
        fc <- (mean(sig_i[[compare]]) - mean(sig_j[[compare]])) /
          mean(sig_j[[compare]])
        conc <- unique(sig_i$Concentration)
      } else {
        fc <- (mean(sig_j[[compare]]) - mean(sig_i[[compare]])) /
                    mean(sig_i[[compare]])
        conc <- unique(sig_j$Concentration)
      }
        
      sig_ij <- wilcox.test(sig_i[[compare]], sig_j[[compare]],
                            alternative = "two.sided")
    
      this_sig <- c(lst[i], lst[j], sig_ij$p.value, fc, conc)
      sig <- rbind(sig, this_sig)
    }
  }
  colnames(sig) <- c("Group1", "Group2", "Significance",
                     "Fold_change", "Concentration")
  rownames(sig) <- c()
  sig <- as.data.frame(sig)
  
  sig$Passage1 <- sapply((strsplit(sig$Group1, "_")), "[[", 1)
  sig$Passage2 <- sapply((strsplit(sig$Group2, "_")), "[[", 1)
  sig <- sig[sig$Passage1 == sig$Passage2, ]
  
  sig$FDR <- p.adjust(sig$Significance, method = "fdr")
  sig$FDR_sig <- ifelse(sig$FDR < 0.05, "Sig", "Non-sig")
  
  sig$Compare <- compare
  return(sig)
}
