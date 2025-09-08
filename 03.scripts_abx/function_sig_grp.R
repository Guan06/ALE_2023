###############################################################################
#  COMPACT LETTER DISPLAY  (base-R only)      ──  with fixed letter for
#                                              a chosen reference population
###############################################################################

make_letters <- function(max_n = 702) {          # a,b,…,z, aa,ab,… (<= 702)
  base <- letters
  if (max_n <= 26) return(base[1:max_n])
  two  <- outer(base, base, paste0)
  c(base, two)[1:max_n]
}

cld <- function(df, alpha = 0.05,
                ref_name = "Parental",          # ← force this pop to “a”
                fixed_letter = "a") {
  
  ## ---- 1. Clean names ------------------------------------------------------
  df$Pop_1 <- trimws(as.character(df$Pop_1))
  df$Pop_2 <- trimws(as.character(df$Pop_2))
  
  pops <- unique(c(df$Pop_1, df$Pop_2))
  
  # If reference exists, place it first so we treat it right away
  if (ref_name %in% pops) {
    pops <- c(ref_name, setdiff(pops, ref_name))
  }
  
  n_pop <- length(pops)
  
  ## ---- 2. Symmetric matrix of p-values ------------------------------------
  pmat  <- matrix(1, n_pop, n_pop, dimnames = list(pops, pops))
  i1    <- match(df$Pop_1, pops)
  i2    <- match(df$Pop_2, pops)
  pmat[cbind(i1, i2)] <- df$P_adj
  pmat[cbind(i2, i1)] <- df$P_adj
  diag(pmat) <- 1
  
  ## ---- 3. Prepare letter matrix -------------------------------------------
  letters_pool <- make_letters()
  
  if (ref_name %in% pops) {
    # Start with one column “a” already assigned to the reference
    letter_mat <- matrix(FALSE, n_pop, 1,
                         dimnames = list(pops, fixed_letter))
    letter_mat[ref_name, fixed_letter] <- TRUE
    next_letter_idx <- 2     # next free label in letters_pool
    # Propagate "a" to every pop not sig. different from the reference
    for (j in seq_len(n_pop)) {
      if (!letter_mat[j, fixed_letter] &&
          pmat[pops[j], ref_name] >= alpha) {
        letter_mat[j, fixed_letter] <- TRUE
      }
    }
  } else {
    letter_mat      <- matrix(FALSE, n_pop, 0, dimnames = list(pops, NULL))
    next_letter_idx <- 1
  }
  
  ## ---- 4. Main greedy assignment with full propagation ---------------------
  for (i in seq_len(n_pop)) {
    pop_i <- pops[i]
    
    ## 4a – try each existing letter
    for (L in seq_len(ncol(letter_mat))) {
      carriers <- pops[letter_mat[, L]]
      if (!letter_mat[pop_i, L] &&              # already has it? skip
          all(pmat[pop_i, carriers] >= alpha))  # compatible?
        letter_mat[pop_i, L] <- TRUE
    }
    
    ## 4b – if still no letter, open a new one
    if (!any(letter_mat[pop_i, ])) {
      letter_mat <- cbind(letter_mat, FALSE)
      colnames(letter_mat)[ncol(letter_mat)] <- letters_pool[next_letter_idx]
      letter_mat[pop_i, ncol(letter_mat)] <- TRUE
      next_letter_idx <- next_letter_idx + 1
    }
    
    ## 4c – PROPAGATE every letter to all further compatible pops
    changed <- TRUE
    while (changed) {
      changed <- FALSE
      for (L in seq_len(ncol(letter_mat))) {
        carriers <- pops[letter_mat[, L]]
        for (j in seq_len(n_pop)) {
          if (!letter_mat[j, L] &&
              all(pmat[pops[j], carriers] >= alpha)) {
            letter_mat[j, L] <- TRUE
            changed <- TRUE
          }
        }
      }
    }
  }
  
  ## ---- 5. Convert logical matrix to compact strings ------------------------
  Group <- apply(letter_mat, 1, function(x)
    paste0(colnames(letter_mat)[x], collapse = ""))
  
  out <- data.frame(Population = pops,
                    Group      = Group,
                    row.names  = NULL,
                    stringsAsFactors = FALSE)
  out[order(out$Population), ]
}

###############################################################################
# Example demonstrating fixed "Parental" letter
###############################################################################
# df_ex <- data.frame(
#   Pop_1 = c("Parental","Parental","Parental","A","A","B"),
#   Pop_2 = c("A","B","C","B","C","C"),
#   P_adj = c(0.30, 0.01, 0.20, 0.04, 0.25, 0.07),  # 0.01 & 0.04 are sig.
#   Antibiotic = NA,
#   Sig        = NA,
#   stringsAsFactors = FALSE
# )


############################################################################
##  compact-letter display with “parental” fixed to letter “a”
############################################################################
assignLetters <- function(dat,
                          alpha          = 0.05,
                          parental_name  = "Parental")
{
  ## ----------------------------------------------------------------------
  ## 1. read the pair-wise table
  ## ----------------------------------------------------------------------
  pops <- unique(c(dat$Pop_1, dat$Pop_2))
  if (! parental_name %in% pops)
    stop(sprintf("Population “%s” not found in the file.", parental_name))
  
  pops <- sort(pops)                      
  k    <- length(pops)
  
  ## ----------------------------------------------------------------------
  ## 2. build the logical matrix  SAME[i,j]  =  TRUE  ⇔  NOT significant
  ## ----------------------------------------------------------------------
  
  SAME <- matrix(TRUE, k, k, dimnames = list(pops, pops))
  for (i in seq_len(nrow(dat))) {
    pos_a <- match(dat$Pop_1[i], pops)      # row index
    pos_b <- match(dat$Pop_2[i], pops)      # column index
    SAME[pos_a, pos_b] <- SAME[pos_b, pos_a] <- (dat$P_adj[i] >= alpha)
  }
  diag(SAME) <- TRUE
  
  ## ----------------------------------------------------------------------
  ## 3. initialise the letter assignments
  ## ----------------------------------------------------------------------
  alphabet <- c(letters, paste0(rep(letters, each = 26), letters))  # a..z, aa, ab, …
  letters_vec  <- setNames(rep("", k), pops)     # storage pop -> letters
  letterGroups <- list()                         # storage letter -> character vector
  
  ## ---- force "parental" to have the letter “a” --------------------------
  letters_vec[parental_name] <- "a"
  letterGroups[["a"]]        <- parental_name
  nextLetter <- 2L           # next free letter index in ‘alphabet’ (1 is “a”)
  
  ## ----------------------------------------------------------------------
  ## 4. helper functions
  ## ----------------------------------------------------------------------
  can_share <- function(pop, group_set)          # TRUE ⇔ pop is non-sig. vs *all* in set
    all(SAME[pop, group_set])
  
  add_letter_to <- function(pop, let) {          # attach letter (if not already there)
    if (!grepl(let, letters_vec[pop]))
      letters_vec[pop] <<- paste0(letters_vec[pop], let)
    letterGroups[[let]] <<- unique(c(letterGroups[[let]], pop))
  }
  
  ## ----------------------------------------------------------------------
  ## 5. main loop through the populations (parental already treated)
  ## ----------------------------------------------------------------------
  for (pop in setdiff(pops, parental_name)) {
    
    assigned <- FALSE
    
    ## (i) try to use an already existing letter
    for (let in names(letterGroups)) {
      if (can_share(pop, letterGroups[[let]])) {
        add_letter_to(pop, let)
        assigned <- TRUE
      }
    }
    
    ## (ii) if nothing fits, create a new letter and
    ##      afterwards give the same new letter to *every*
    ##      earlier population that can share it
    if (!assigned) {
      let <- alphabet[nextLetter];  nextLetter <- nextLetter + 1L
      letterGroups[[let]] <- character(0)        # create the slot
      add_letter_to(pop, let)                    # give new letter to current pop
      
      ## try to add the new letter to previous pops where possible
      for (old in pops[pops != pop]) {
        if (can_share(old, letterGroups[[let]]))
          add_letter_to(old, let)
      }
    }
  }
  
  ## ----------------------------------------------------------------------
  ## 6. sort letters inside each label (a before b …) and return
  ## ----------------------------------------------------------------------
  letters_vec <- vapply(letters_vec,
                        function(x) paste(sort(unique(strsplit(x,"")[[1]])),
                                          collapse = ""),
                        character(1))
  out <- data.frame(Population = names(letters_vec),
                    Letter     = unname(letters_vec),
                    row.names  = NULL,
                    stringsAsFactors = FALSE)
  
  out[order(out$Population), ]
}

###############################################################################
# Example demonstrating fixed "Parental" letter
###############################################################################
##############################
## save this as  pairs.tsv
# Pop_1  Pop_2  P_adj
# parental A     0.690
# parental B     0.010
# parental C     0.002
# A        B     0.005
# A        C     0.730
# B        C     0.004
##########################
# assignLetters("pairs.tsv")
#   Population Letter
# 1          A     ac
# 2          B      b
# 3          C      c
# 4   parental      a
