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







#Analysed the compact model
model_compact = readRDS("results/varComp/model_compact_results.rds") %>% tbl_df() %>%
  dplyr::filter(converged == TRUE) %>%
  dplyr::select(-type, -converged) %>% #Remove unnecessary columns
  dplyr::rename(residual = Residual,SL1344_IFNg = `SL1344:IFNg`, batch = salmonella)

#Bin by residual and calculate variance explained within each bin
binned_table = binGenesByResidual(model_compact, n_bins = 20)
var_explained = meanVarianceWithinBins(binned_table)
var_comp_plot = plotBinnedVariance(var_explained) 
ggsave("results/varComp/var_comp_plot.pdf", plot = var_comp_plot, width = 11, height = 7)

#Make a violin plot of all components
dat = tidyr::gather(model_compact, factor, var_explained, donor:SL1344_IFNg)
ggplot(dat, aes(x = factor, y = var_explained)) + geom_violin(scale = "width" ) + geom_boxplot(width = .2, outlier.shape = NA)

#### Pool technical source of variance together ####
model_pooled = dplyr::transmute(model_compact, gene_id, donor, gender, IFNg, SL1344, SL1344_IFNg, residual, 
                                technical = batch + library_pool + rna_concentration, purity = purity_bins)
binned_pooled = binGenesByResidual(model_pooled, n_bins = 20)
var_pooled = meanVarianceWithinBins(binned_pooled)
var_pooled_plot = plotBinnedVariance(var_pooled)
ggsave("results/varComp/var_pooled_plot.pdf", plot = var_pooled_plot, width = 11, height = 7)

#Bin variance by maximum factor
maximum_factor = maximumFactorPerGene(model_compact)
ggplot(maximum_factor, aes(x = component_max, y = var_explained_max)) + geom_violin()
table(maximum_factor$component_max)

#Calculate the number of genes in each bin
bin_size_table = dplyr::group_by(binned_table, residual_bin) %>% dplyr::summarise(bin_size = length(residual_bin))

#Make plots
bin_sizes_plot = ggplot(bin_size_table, aes(x = residual_bin, y = bin_size)) + geom_bar(stat = "identity") + 
  ylab("Number of genes") + 
  xlab("Residual variance bin")
ggsave("results/varComp/bin_sizes_plot.pdf", plot = bin_sizes_plot, width = 8, height = 3)

