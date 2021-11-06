### Gene Expression Data ###

# Data from https://github.com/genomicsclass/tissuesGeneExpression
devtools::install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)
data(tissuesGeneExpression)

data_gene <- t(e)
data_gene <- as.data.frame(data_gene)
