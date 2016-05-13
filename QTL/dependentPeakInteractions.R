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
naive_sl1344_atac = extractConditionFromExpressionList(c("naive","SL1344"), atac_data)
naive_ifng_sl1344_atac = extractConditionFromExpressionList(c("naive","IFNg_SL1344"), atac_data)

peak_pairs = data_frame(master_id = "ATAC_peak_145162", dependent_id = "ATAC_peak_145165", snp_id = "rs10928660")

#Specify models of interest
model1 = as.formula("cqn ~genotype + peak_type + condition_name + cqn_PC1 + cqn_PC2 + cqn_PC3 + peak_type*condition_name + 
               genotype*peak_type + genotype*condition_name + genotype*condition_name*peak_type")
model0 = as.formula("cqn ~genotype + peak_type + condition_name  + cqn_PC1 + cqn_PC2 + cqn_PC3 + peak_type*condition_name + 
                    genotype*peak_type + genotype*condition_name")

#Test for threeway interaction
res = testThreewayInteraction(peak_pairs, naive_ifng_atac$cqn, naive_ifng_atac$sample_metadata, vcf_file, model0, model1, p_only = FALSE)
ggplot(res$data, aes(x = factor(genotype), y = cqn)) + geom_point() + facet_grid(peak_type ~ condition_name)

#Import list of dependent peaks
dependent_peaks = readRDS("results/ATAC/QTLs/depdendent_peaks.rds")

#Test peak-peak-genotype interactions in naive_vs_IFNg setting
peak_peak_interactions = purrr::by_row(dependent_peaks$unique_masters, testThreewayInteraction, naive_ifng_atac$cqn, 
                  naive_ifng_atac$sample_metadata, vcf_file, model0, model1, .collate = "rows", .to = "naive_vs_IFNg_pvalue")
peak_peak_inter = dplyr::mutate(peak_peak_interactions, naive_vs_IFNg_fdr = p.adjust(naive_vs_IFNg_pvalue, method = "fdr")) %>% 
  dplyr::arrange(naive_vs_IFNg_fdr) %>% 
  dplyr::filter(naive_vs_IFNg_fdr < 0.1)

peak_peak_interactions_sl = purrr::by_row(dependent_peaks$unique_masters, testThreewayInteraction, naive_sl1344_atac$cqn, 
                                       naive_sl1344_atac$sample_metadata, vcf_file, model0, model1, .collate = "rows", .to = "naive_vs_SL1344_pvalue")
peak_peak_inter_sl = dplyr::mutate(peak_peak_interactions_sl, naive_vs_SL1344_fdr = p.adjust(naive_vs_SL1344_pvalue, method = "fdr")) %>% 
  dplyr::arrange(naive_vs_SL1344_fdr) %>% 
  dplyr::filter(naive_vs_SL1344_fdr < 0.1)

#Visualise a couple of examples
res = testThreewayInteraction(peak_peak_inter[3,], naive_ifng_atac$cqn, naive_ifng_atac$sample_metadata, vcf_file, model0, model1, p_only = FALSE)
ggplot(res$data, aes(x = factor(genotype), y = cqn)) + geom_point() + geom_boxplot() + facet_grid(peak_type ~ condition_name)

#Make plots
data = purrr::by_row(peak_peak_inter, testThreewayInteraction, naive_ifng_atac$cqn, 
              naive_ifng_atac$sample_metadata, vcf_file, model0, model1, p_only = FALSE)
data_list = map(as.list(data$.out), function(x){x$data})
plot_list = map(data_list, ~ggplot(., aes(x = factor(genotype), y = cqn)) + geom_point() + geom_boxplot() + facet_grid(peak_type ~ condition_name))
names(plot_list) = paste(data$dependent_id, data$master_id, sep = "-")

data = purrr::by_row(peak_peak_inter_sl, testThreewayInteraction, naive_sl1344_atac$cqn, 
                     naive_sl1344_atac$sample_metadata, vcf_file, model0, model1, p_only = FALSE)
data_list = map(as.list(data$.out), function(x){x$data})
plot_list = map(data_list, ~ggplot(., aes(x = factor(genotype), y = cqn)) + geom_point() + geom_boxplot() + facet_grid(peak_type ~ condition_name))
names(plot_list) = paste(data$dependent_id, data$master_id, sep = "-")

savePlotList(plot_list, "results/ATAC/QTLs/peak-peak_interactions/naive_vs_SL1344/", suffix = ".pdf", width = 8, height = 8)


##### Distances between master and dependent peaks ####
#Calculated distances between master and dependent peaks
peak_distances = calculatePeakDistance(dependent_peaks$unique_masters, atac_data$gene_metadata) %>%
  dplyr::select(dependent_id, master_id, distance)
peak_distances_cl = calculatePeakDistance(dependent_peaks$cluster_master, atac_data$gene_metadata) %>%
  dplyr::select(dependent_id, master_id, distance)

#Make a histogram of distances
peak_distance_df = rbind(peak_distances, peak_distances_cl)
plot = ggplot(peak_distance_df, aes(x = distance/1000)) + geom_histogram(binwidth = 1) + 
  xlab("Distance between master and dependent peaks (kb)") + 
  ylab("Number of peak pairs")
ggsave("results/ATAC/QTLs/properties/master_dependent_peak_distance.pdf", plot = plot, width = 6, height = 6)


#Calculate distances within peak clusters


