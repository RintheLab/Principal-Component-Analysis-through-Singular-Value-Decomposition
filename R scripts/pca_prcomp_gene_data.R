### PCA using prcomp() function on expression data ###

library(ggplot2)
library(purrr)
library(dplyr)
library(ggpubr)

# Data --------------------------------------------------------------------

source("R scripts/gene_expression_data.R")

# PCA with prcomp ---------------------------------------------------------

dg_pca <- prcomp(data_gene)

# Scree plot --------------------------------------------------------------

# Eigenvalues
ev_dg_prcomp <- dg_pca$sdev^2 / (nrow(data_gene) - 1)

# Percentage of variation for each PC
per_var_prcomp <- data.frame(ev_dg_prcomp) %>% 
  mutate(
    per_var = ev_dg_prcomp * 100 / sum(ev_dg_prcomp),
    pc = map_chr(
      1:nrow(dg_pca$x), .f = function(x) paste0("PC", as.character(x))
      )
  )

# Scree plot with the first 15 principal components 
scree_plot_dg_prcomp <- ggplot(
  per_var_prcomp[1:15,], 
  aes(reorder(pc, -per_var), per_var)
) +
  geom_col() +
  xlab("Principal component") +
  ylab("Percentage of variation (%)") +
  theme_classic()


# Score plot --------------------------------------------------------------

# Scores
scores_dg_prcomp <- dg_pca$x %>% 
  as.data.frame() %>% 
  mutate(Tissue = tissue) %>% 
  relocate(Tissue)

# Score plot
score_plot_prcomp <- ggplot(scores_dg_prcomp, aes(PC1, PC2, color = Tissue)) +
  geom_point(size = 2) +
  xlab("PC1 (33%)") + 
  ylab("PC2 (14%)") +
  ggtitle("PCA on expression data using prcomp()") +
  theme_classic()


# Loadings ----------------------------------------------------------------

# Calculate loadings
ls_dg_prcomp <- dg_pca$rotation %*% diag(dg_pca$sdev)

dim(ls_dg_prcomp)

# First 5 rows and columns
ls_dg_prcomp[1:5, 1:5]


# Eigenvectors (usually reported as loadings) -----------------------------

dg_pca$rotation[1:5, 1:5]
