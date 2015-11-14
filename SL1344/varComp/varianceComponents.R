#This script performs the timecomsuming bit of fitting mixed-effects
#models on normalized gene expression data
library("plyr")
library("dplyr")
library("lme4")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")

#Import eQTL data list
data_list = readRDS("results/SL1344/eqtl_data_list.rds")

#Perform PCA
pca_list = performPCA(data_list$exprs_cqn, data_list$sample_metadata)
pca_plot = ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, 
                                           color = condition_name, label = sample_id)) + 
  geom_text()

#Make some modifications to the sample_metadata (TODO: move this code to line_metadata construction)
data_list$sample_metadata = dplyr::mutate(data_list$sample_metadata,
                                          rna_auto = ifelse(chemistry == "V4_auto", "yes", "no"))

#Prepare data for lmer analysis
gene_ids = rownames(data_list$exprs_cqn)
gene_id_list = idVectorToList(gene_ids)
gene_data_list = lapply(gene_id_list, constructGeneData, 
                        data_list$exprs_cqn, data_list$sample_metadata)

#Fit the extended model
model_extended <- function(model_data){
  mod = lme4::lmer(value ~ (1|SL1344) + (1|IFNg) + (1|SL1344:IFNg) + 
                     (1|salmonella) + (1|donor) + (1|gender) + (1|purity_bins) +
                     (1|rna_concentration) + (1|passage_diff_bins) + (1|diff_bins) + 
                     (1|library_pool) + (1|mf_diff_days) + (1|rna_auto), model_data)
  return(mod)
}

variance_list = lapply(gene_data_list, estimateVarianceExplained, model_extended)
var_table = ldply(variance_list, .id = "gene_id")
saveRDS(var_table, "results/SL1344/varComp/model_extended_results.rds")

#Fit the compact model
#model_compact <- function(model_data){
#  mod = lme4::lmer(value ~ (1|SL1344) + (1|IFNg) + (1|SL1344:IFNg) + 
#                     (1|salmonella) + (1|donor) + (1|gender) + (1|purity_bins) +
#                     (1|rna_concentration) + (1|library_pool) + (1|rna_auto), model_data)
#  return(mod)
#}

#variance_list = lapply(gene_data_list, estimateVarianceExplained, model_compact)
#var_table = ldply(variance_list, .id = "gene_id")
#saveRDS(var_table, "results/SL1344/varComp/model_compact_results.rds")

