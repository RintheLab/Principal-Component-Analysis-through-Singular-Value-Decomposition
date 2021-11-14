### PCA using svd() function on expression data ###

library(ggplot2)
library(purrr)
library(dplyr)
library(ggpubr)

# Data --------------------------------------------------------------------

source("R scripts/gene_expression_data.R")

head(data_gene)[, 1:5]

# Centered data

gene_means <- map_dbl(data_gene, mean)

data_gene_c <- map2(data_gene, gene_means, .f = function(x, mean) x - mean)
data_gene_c <- data.frame(data_gene_c)

# Singular value decomposition --------------------------------------------

dg_svd <- svd(data_gene_c)

# U matrix 
dg_u <- dg_svd$u

# V matrix
dg_v <- dg_svd$v

# Singular values
dg_d <- diag(dg_svd$d)

# Scree plot --------------------------------------------------------------

# Eigenvalues
ev_dg_svd <- dg_svd$d^2 / (nrow(data_gene) - 1)

# Percentage of variation for each PC
pervar_dg_svd <- data.frame(ev_dg_svd) %>% 
  mutate(
    per_var = ev_dg_svd * 100 / sum(ev_dg_svd),
    pc = map_chr(1:nrow(dg_u), .f = function(x) paste0("PC", as.character(x))
    )
  )

# Scree plot with the first 15 principal components 
scree_plot_dg_svd <- ggplot(
  pervar_dg_svd[1:15,], 
  aes(reorder(pc, -per_var), per_var)
) +
  geom_col() +
  xlab("Principal component") +
  ylab("Percentage of variation (%)") +
  theme_classic()


# Score plot --------------------------------------------------------------

# Scores 
scores_dg_svd <- data.frame(dg_u %*% dg_d)

# Change the default names for PC names
colnames(scores_dg_svd) <- map_chr(
  1:nrow(scores_dg_svd), .f = function(x) paste0("PC", as.character(x))
)

# Add tissue column
scores_dg_svd <- scores_dg_svd %>% 
  mutate(Tissue = tissue) %>% 
  relocate(Tissue)

# Score plot
score_plot_svd <- ggplot(scores_dg_svd, aes(PC1, PC2, color = Tissue)) +
  geom_point(size = 2) +
  xlab("PC1 (33%)") + 
  ylab("PC2 (14%)") +
  ggtitle("PCA on gene expression data using svd()") +
  theme_classic()


# Loadings ----------------------------------------------------------------

# Calculate loadings
ls_dg_svd <- (dg_v %*% dg_d) / sqrt(nrow(data_gene) - 1)

dim(ls_dg_svd)

# First 5 rows and columns
ls_dg_svd[1:5, 1:5]


# Eigenvectors (usually reported as loadigins) ----------------------------

dg_v[1:5, 1:5]
