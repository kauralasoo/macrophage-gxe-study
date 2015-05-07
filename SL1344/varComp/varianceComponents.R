#This script performs the timecomsuming bit of fitting mixed-effects models on normalized gene expression data
library("plyr")
library("dplyr")
library("lme4")
library("devtools")
library("ggplot2")
load_all("macrophage-gxe-study/seqUtils/")

#Import line metadata
line_metadata = readRDS("results/SL1344/sample_info/compiled_line_metadata.rds")
expression_list = readRDS("results/SL1344/combined_expression_data.rds")

#Remove fpdj_2 sample
design = dplyr::filter(expression_list$design, !(donor == "fpdj" & replicate == 2))
rownames(design) = design$sample_id

#Calculate mean expression per condition
exprs_cqn = expression_list$exprs_cqn[,rownames(design)]
exprs_cqn_mean = calculateMean(expression_list$exprs_cqn, expression_list$design, "condition")
expressed_genes = names(which(apply(exprs_cqn_mean, 1, max) > 1))
exprs_cqn_filtered = exprs_cqn[expressed_genes,]

#Perform PCA
pca_list = performPCA(exprs_cqn_filtered, design)
pca_plot = ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition_name, label = sample_id)) + 
  geom_text()

#Specify alternative models
model_extended <- function(model_data){
  mod = lme4::lmer(exp_value ~ (1|SL1344) + (1|IFNg) + (1|SL1344:IFNg) + 
                     (1|salmonella) + (1|donor) + (1|gender) + (1|purity_bins) +
                     (1|rna_concentration) + (1|passage_diff_bins) + (1|diff_bins) + 
                     (1|library_pool) + (1|mf_diff_days), model_data)
  return(mod)
}

model_compact <- function(model_data){
  mod = lme4::lmer(exp_value ~ (1|SL1344) + (1|IFNg) + (1|SL1344:IFNg) + 
                     (1|salmonella) + (1|donor) + (1|gender) + (1|purity_bins) +
                     (1|rna_concentration) + (1|library_pool), model_data)
  return(mod)
}

#Apply lmer to a list of genes
gene_ids = rownames(exprs_cqn_filtered)
gene_id_list = idVectorToList(gene_ids)
gene_data_list = lapply(gene_id_list, constructGeneData, exprs_cqn, design, line_metadata)

#Fit the compact model
variance_list = lapply(gene_data_list, estimateVarianceExplained, model_compact)
var_table = ldply(variance_list, .id = "gene_id")
saveRDS(var_table, "results/varComp/model_compact_results.rds")

#Fit the extended model
variance_list = lapply(gene_data_list, estimateVarianceExplained, model_extended)
var_table = ldply(variance_list, .id = "gene_id")
saveRDS(var_table, "results/varComp/model_extended_results.rds")
