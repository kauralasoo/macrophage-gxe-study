library("purrr")
library("ggplot2")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Import ATAC data
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import minimal p-values
min_pvalue_list = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(min_pvalue_list, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name")

#Extract data
naive_ifng_atac = extractConditionFromExpressionList(c("naive","IFNg"), atac_data)
peak_pairs = data_frame(master_id = "ATAC_peak_145162", dependent_id = "ATAC_peak_145165", snp_id = "rs10928660")

#Specify models of interest
model1 = as.formula("cqn ~genotype + peak_type + condition_name + cqn_PC1 + cqn_PC2 + cqn_PC3 + peak_type*condition_name + 
               genotype*peak_type + genotype*condition_name + genotype*condition_name*peak_type")
model0 = as.formula("cqn ~genotype + peak_type + condition_name  + cqn_PC1 + cqn_PC2 + cqn_PC3 + peak_type*condition_name + 
                    genotype*peak_type + genotype*condition_name")

#Test for threeway interaction
res = testThreewayInteraction(peak_pairs, naive_ifng_atac$cqn, naive_ifng_atac$sample_metadata, vcf_file, model0, model1, p_only = FALSE)
ggplot(res$data, aes(x = factor(genotype), y = cqn)) + geom_point() + facet_grid(peak_type ~ condition_name)

