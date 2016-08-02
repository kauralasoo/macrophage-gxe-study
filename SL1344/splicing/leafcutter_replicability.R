library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("Rsamtools")
library("purrr")

#Import data
prop_list = readRDS("results/SL1344/combined_proportions.row_quantile.rds")

#Import variant information
snp_info = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

#Find minimal p-values from fastQTL results
naive_fqtl = importFastQTLTable("results/SL1344/leafcutter/fastqtl_output/naive_100kb_permuted.txt.gz")
ifng_fqtl = importFastQTLTable("results/SL1344/leafcutter/fastqtl_output/IFNg_100kb_permuted.txt.gz")
sl1344_fqtl = importFastQTLTable("results/SL1344/leafcutter/fastqtl_output/SL1344_100kb_permuted.txt.gz")
ifng_sl1344_fqtl = importFastQTLTable("results/SL1344/leafcutter/fastqtl_output/IFNg_SL1344_100kb_permuted.txt.gz")

fastqtl_pvalue_list = list(naive = naive_fqtl,
                           IFNg = ifng_fqtl,
                           SL1344 = sl1344_fqtl, 
                           IFNg_SL1344 = ifng_sl1344_fqtl)
saveRDS(fastqtl_pvalue_list, "results/SL1344/leafcutter/leafcutter_min_pvalues.rds")

#Calculate Pi1
pi1_stat = calculatePairwisePi1(fastqtl_pvalue_list)
write.table(pi1_stat, "results/SL1344/leafcutter/fastqtl_pi1_results.txt", sep = "\t", quote = FALSE)
pi1_stat_tidy = calculatePairwisePi1(fastqtl_pvalue_list, tidy = TRUE)
write.table(pi1_stat_tidy, "results/SL1344/leafcutter/fastqtl_pi1_results_tidy.txt", sep = "\t", quote = FALSE)

#Choose one intron randomly per cluster
cluster_ids = dplyr::select(prop_list$gene_metadata, gene_id, cluster_id)
randon_intron = purrr::map(fastqtl_pvalue_list, ~dplyr::left_join(., cluster_ids, by = "gene_id") %>% 
                 dplyr::arrange(gene_id) %>%
                 dplyr::group_by(cluster_id) %>%
                 dplyr::filter(row_number() == 1))

#Calculate PI1
pi1_stat = calculatePairwisePi1(randon_intron)
write.table(pi1_stat, "results/SL1344/leafcutter/fastqtl_rnd_pi1_results.txt", sep = "\t", quote = FALSE)
pi1_stat_tidy = calculatePairwisePi1(randon_intron, tidy = TRUE)
write.table(pi1_stat_tidy, "results/SL1344/leafcutter/fastqtl_rnd_pi1_results_tidy.txt", sep = "\t", quote = FALSE)


 
