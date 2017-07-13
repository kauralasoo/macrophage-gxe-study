library("SummarizedExperiment")
library("ggplot2")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Import Expression data
se_fairfax = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")

#Extract data
sample_metadata = tbl_df2(colData(se_fairfax))
exprs_matrix = assays(se_fairfax)$exprs

#calculate mean expression
mean = calculateMean(exprs_matrix, as.data.frame(sample_metadata), "condition_name")

#Use PCA to id outliers
pca_res = performPCA(exprs_matrix, as.data.frame(sample_metadata), n_pcs = 6)

#Make some PC plots
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = condition_name)) + geom_point()
ggplot(pca_res$pca_matrix, aes(x = PC2, y = PC3, color = condition_name)) + geom_point()