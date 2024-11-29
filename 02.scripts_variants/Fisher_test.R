## For all detected variants, AF > 0.05
# Number of samples and variants in ALE and ALE 2.0
data <- matrix(c(12, 323, 286, 10102), nrow = 2)

# Name the rows and columns (optional, for clarity)
rownames(data) <- c("ALE_2.0","ALE")
colnames(data) <- c("Sample","Variant")

# Perform Fisher's exact test
fisher.test(data)

## For predominant variants, AF > 0.5
# Number of samples and variants in ALE and ALE 2.0
data <- matrix(c(12, 323, 26, 310), nrow = 2)

# Name the rows and columns (optional, for clarity)
rownames(data) <- c("ALE_2.0","ALE")
colnames(data) <- c("Sample","Variant")

# Perform Fisher's exact test
fisher.test(data)
