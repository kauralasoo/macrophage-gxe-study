library("readr")
load_all("../seqUtils/")

#Load RASQUAL p-values
rasqual_pvalues = seqUtils::importRasqualTable("results/ATAC/rasqual/output/naive_50kb.filtered.txt")
rasqual_pvalues = seqUtils::importRasqualTable("results/ATAC/rasqual/output/naive_50kb_ls.txt.gz")

rasqual_min_pvalues = dplyr::group_by(rasqual_pvalues, gene_id) %>%
  dplyr::arrange(p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::mutate(p_bonferroni = p.adjust(p_nominal, "bonferroni", n_cis_snps)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(p_nominal) 

rasqual_bonferroni_filtered = dplyr::filter(rasqual_min_pvalues, p_bonferroni < 0.05)

#Load FastQTL p-values
naive_perm_pvals = importFastQTLTable("results/ATAC/fastqtl/output/naive_50kb_tpm_permuted.txt.gz") %>% tbl_df()
p_vals = dplyr::group_by(naive_perm_pvals, gene_id) %>% 
  dplyr::mutate(p_bonferroni = p.adjust(p_nominal, "bonferroni", n_cis_snps)) %>%
  dplyr::filter(gene_id %in% rasqual_min_pvalues$gene_id) %>%
  dplyr::filter(snp_id %in% unique(rasqual_pvalues$snp_id)) %>%
  ungroup() %>%
  dplyr::mutate(p_fdr = p.adjust(p_beta, "fdr"))

fastqtl_bonferroni_filtered = dplyr::filter(p_vals, p_bonferroni < 0.05) %>% 
  dplyr::arrange(p_bonferroni)

fastqtl_fdr_filtered = dplyr::filter(p_vals, p_fdr < 0.05) %>% 
  dplyr::arrange(p_fdr)

hist(dplyr::filter(rasqual_min_pvalues, gene_id %in% fastqtl_bonferroni_filtered$gene_id)$p_nominal,  breaks  = 50)
hist(dplyr::filter(p_vals, gene_id %in% rasqual_bonferroni_filtered$gene_id)$p_nominal, breaks  = 50)

list_union = union(rasqual_bonferroni_filtered$gene_id, fastqtl_bonferroni_filtered$gene_id)
snp_union = union(rasqual_bonferroni_filtered$snp_id, fastqtl_bonferroni_filtered$snp_id)

r = dplyr::filter(rasqual_pvalues, gene_id %in% list_union)
f = dplyr::filter(naive_perm_pvals, gene_id %in% list_union, snp_id %in% snp_union)
d = dplyr::left_join(r, f, by = c("gene_id","snp_id"))
ggplot2::ggplot(d, aes(x = -log10(p_nominal.x), y=-log10(p_nominal.y) )) +
  geom_point() +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0, 20))

