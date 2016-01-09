library("readr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")

#Import ATAC data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
peak_positions = dplyr::select(atac_list$gene_metadata, gene_id, start, end) %>%
  dplyr::mutate(centre = start + floor((end - start)/2))

#Load RASQUAL p-values
rasqual_pvalues = seqUtils::importRasqualTable("results/ATAC/rasqual/output/naive_chr11_50kb.txt")
rasqual_pvalues_cov = seqUtils::importRasqualTable("results/ATAC/rasqual/output/naive_chr11_50kb_cov.txt.gz")
rasqual_pvalues_cov_pc4 = seqUtils::importRasqualTable("results/ATAC/rasqual/output/naive_chr11_50kb_cov_pc4.txt.gz")
rasqual_pvalues_gc = seqUtils::importRasqualTable("results/ATAC/rasqual/output/naive_chr11_50kb_gc.txt.gz")

rasqual_min_pvalues = findMinimalSnpPvalues(rasqual_pvalues)
rasqual_min_pvalues_cov = findMinimalSnpPvalues(rasqual_pvalues_cov)
rasqual_min_pvalues_cov_pc4 = findMinimalSnpPvalues(rasqual_pvalues_cov_pc4)
rasqual_min_pvalues_gc = findMinimalSnpPvalues(rasqual_pvalues_gc)

rasqual_fdr_filtered = dplyr::filter(rasqual_min_pvalues, p_fdr < 0.1)
dplyr::filter(rasqual_min_pvalues_cov_pc4, p_fdr < 0.1)

#Load FastQTL p-values
naive_perm_pvals = importFastQTLTable("results/ATAC/fastqtl/output/naive_chr11_50kb_tpm_perm.txt.gz") %>% tbl_df()
fastqtl_min_pvalues = findMinimalSnpPvalues(naive_perm_pvals) %>%
  dplyr::mutate(p_fdr_beta = p.adjust(p_beta, "fdr"))

fastqtl_fdr_filtered = dplyr::filter(fastqtl_min_pvalues, p_fdr < 0.1)
dplyr::filter(fastqtl_min_pvalues, p_fdr_beta < 0.05)

hist(dplyr::filter(rasqual_min_pvalues, gene_id %in% fastqtl_bonferroni_filtered$gene_id)$p_nominal,  breaks  = 50)
hist(dplyr::filter(p_vals, gene_id %in% rasqual_bonferroni_filtered$gene_id)$p_nominal, breaks  = 50)

list_union = union(rasqual_fdr_filtered$gene_id, fastqtl_fdr_filtered$gene_id)
snp_union = union(rasqual_bonferroni_filtered$snp_id, fastqtl_bonferroni_filtered$snp_id)

r = dplyr::filter(rasqual_min_pvalues, gene_id %in% list_union)
f = dplyr::filter(naive_perm_pvals, gene_id %in% list_union)
d = dplyr::left_join(r, f, by = c("gene_id"))
ggplot2::ggplot(d, aes(x = -log10(p_nominal.x), y=-log10(p_nominal.y) )) +
  geom_point() +
  scale_x_continuous(limits = c(0, 25)) +
  scale_y_continuous(limits = c(0, 25))

#Compare 50kb vs 500kb cis windows
rasqual_pvalues = seqUtils::importRasqualTable("results/ATAC/rasqual/output/naive_chr11_50kb.txt.gz")
rasqual_min_pvalues = findMinimalSnpPvalues(rasqual_pvalues)
rasqual_50kb_hits = dplyr::left_join(rasqual_min_pvalues, peak_positions, by = "gene_id") %>%
  dplyr::mutate(distance = pos - centre)
rasqual_fdr_filtered = dplyr::filter(rasqual_50kb_hits, p_fdr < 0.1)

rasqual_pvalues_500kb = seqUtils::importRasqualTable("results/ATAC/rasqual/output/naive_chr11_500kb.txt.gz")
rasqual_min_pvalues_500kb = findMinimalSnpPvalues(rasqual_pvalues_500kb)
rasqual_500kb_hits = dplyr::left_join(rasqual_min_pvalues_500kb, peak_positions, by = "gene_id") %>%
  dplyr::mutate(distance = pos - centre)
rasqual_fdr_filtered_500kb = dplyr::filter(rasqual_500kb_hits, p_fdr < 0.1)


