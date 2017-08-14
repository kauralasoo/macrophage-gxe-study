library("readr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("ggplot2")
library("GenomicRanges")

#Filter ASB counts
filterASBCounts <- function(ASB_counts, total_count = 10, min_count = 2){
  counts = dplyr::transmute(ASB_counts, chr = contig, pos = position, snp_id = variantID, refAllele, altAllele, refCount, altCount, totalCount) %>%
    dplyr::filter(totalCount >= total_count) %>%
    dplyr::mutate(minCount = pmin(refCount, altCount)) %>%
    dplyr::filter(minCount >= min_count) %>%
    dplyr::mutate(ratio = altCount/totalCount)
  return(counts)
}

#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Import QTL credible sets
credible_sets_df = importCredibleSets("results/ATAC/QTLs/rasqual_credible_sets.rds", atac_list$gene_metadata)
naive_effects = dplyr::filter(credible_sets_df$naive, gene_id == overlap_peak_id) %>%
  dplyr::transmute(peak_id = gene_id, snp_id, beta, effect_size, p_nominal)

#Import counts
ASB_counts = readr::read_tsv("processed/Schultze/ASEcounts/PU1_naive.ASEcounts.gz", col_names = TRUE, col_types = "ciccciiiiiiii")
pu1_counts = filterASBCounts(ASB_counts, total_count = 10, min_count = 2)
  
#Compare caQTL effect to effect on biding
filtered_pu1 = dplyr::left_join(pu1_counts, naive_effects) %>% 
  dplyr::filter(!is.na(beta))

t = cor.test(filtered_pu1$ratio, filtered_pu1$effect_size, method = "spearman")
pu1_plot = ggplot(filtered_pu1, aes(x = ratio, y = effect_size)) + 
  geom_point() + 
  ylab(expression(paste("RASQUAL effect size (",pi,")"))) +
  xlab(expression(paste("PU.1 allelic ratio"))) +
  theme_light() +
  annotate("text",x = 0.25, y = 0.92, label = paste0("rho = ", round(t$estimate,2)))
ggsave("figures/supplementary/PU1_ASB.pdf", plot = pu1_plot, width = 5, height = 5)


#CEBP_beta overlaps
ASB_counts = readr::read_tsv("processed/OCallaghan/ASEcounts/CEBPbeta_ctrl_204.ASEcounts.gz", col_names = TRUE, col_types = "ciccciiiiiiii")
cebpb_counts = filterASBCounts(ASB_counts, total_count = 10, min_count = 2, min_ratio = 0.1)

#Compare caQTL effect to effect on biding
filtered_cebpb = dplyr::left_join(cebpb_counts, naive_effects) %>% 
  dplyr::filter(!is.na(beta))

t = cor.test(filtered_cebpb$ratio, filtered_cebpb$effect_size, method = "spearman")
cebpb_plot = ggplot(filtered_cebpb, aes(x = ratio, y = effect_size)) + 
  geom_point() + 
  ylab(expression(paste("RASQUAL effect size (",pi,")"))) +
  xlab(expression(paste("CEBP", beta, " allelic ratio"))) +
  theme_light() +
  annotate("text",x = 0.25, y = 0.85, label = paste0("rho = ", round(t$estimate,2)))
ggsave("figures/supplementary/CEBPb_ASB.pdf", plot = cebpb_plot, width = 5, height = 5)


#Import caQTL and eQTL pairs
caQTL_eQTL_pairs = readRDS("results/ATAC_RNA_overlaps/caQTL_eQTL_pairs_betas.rds")
filtered_pairs = dplyr::filter(caQTL_eQTL_pairs, condition_name == "naive", phenotype == "ATAC-seq")

#Identify two groups of peaks
appear_peaks = dplyr::filter(filtered_pairs, abs(beta) < 0.59) %>% dplyr::select(peak_id) %>% unique()
persistent_peaks = dplyr::filter(filtered_pairs, abs(beta) > 0.59) %>% dplyr::select(peak_id) %>% unique()


