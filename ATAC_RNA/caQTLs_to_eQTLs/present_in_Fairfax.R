library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")

priming_pairs = readRDS("results/ATAC_RNA_overlaps/caQTL_eQTL_pairs_betas.rds")

priming_peaks = dplyr::filter(priming_pairs, type == "Indirect", phenotype == "RNA-seq") %>% dplyr::select(peak_id, gene_id, snp_id) %>% dplyr::distinct()

#Import eQTL p-values
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_pvalues = dplyr::filter(rasqual_min_pvalues$naive, p_fdr < fdr_thresh) %>% dplyr::transmute(gene_id, eQTL_snp_id = snp_id)

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")


my_eQTLs = dplyr::left_join(priming_peaks, rasqual_pvalues, by = "gene_id") %>% 
  dplyr::filter(!is.na(eQTL_snp_id)) %>%
  dplyr::group_by(peak_id, gene_id, snp_id) %>%
  dplyr::mutate(R2 = calculatePairR2(snp_id, eQTL_snp_id, vcf_file$genotypes)) %>%
  dplyr::filter(R2 > 0.8)

#Visualize these eQTLs
a = dplyr::filter(priming_pairs, gene_id %in% my_eQTLs$gene_id)
plot = ggplot(a, aes(x = stimulation_state, y = beta, group = gene_id, color = gene_name)) + geom_point() + geom_line() + facet_wrap(~phenotype) +
  ylab("Log2 fold change")
ggsave("figures/supplementary/naive_eQTL_enhancer_priming.png", plot = plot, width = 5, height = 4)

  
#Repeat analysis for Fairfax data
fairfax_min_pvalues = readRDS("results/Fairfax/fairfax_qtl_min_pvalues.rds")
fairfax_pvalues = dplyr::filter(fairfax_min_pvalues$full$CD14, p_fdr < 0.1) %>% dplyr::transmute(gene_id = group_id, eQTL_snp_id = snp_id)

fairfax_QTLs = dplyr::left_join(priming_peaks, fairfax_pvalues, by = "gene_id") %>% 
  dplyr::filter(!is.na(eQTL_snp_id)) %>%
  dplyr::filter(eQTL_snp_id %in% vcf_file$snpspos$snpid) %>%
  dplyr::group_by(peak_id, gene_id, snp_id) %>%
  dplyr::mutate(R2 = calculatePairR2(snp_id, eQTL_snp_id, vcf_file$genotypes)) %>%
  dplyr::filter(R2 > 0.8)


#Look at the direction of effect size in the stimulated condition
priming_pairs = readRDS("results/ATAC_RNA_overlaps/caQTL_eQTL_pairs_betas.rds")
sign_count = dplyr::filter(priming_pairs, stimulation_state == "Stimulated") %>% 
  dplyr::select(gene_id, peak_id, snp_id, phenotype, beta, max_effect) %>% 
  dplyr::mutate(phenotype = ifelse(phenotype == "ATAC-seq", "ATAC", "RNA")) %>% 
  dplyr::distinct() %>% tidyr::spread(phenotype, beta) %>%
  dplyr::mutate(diff_sign = sign(ATAC)-sign(RNA))


#When using a more stringent threshold
priming_pairs = readRDS("results/ATAC_RNA_overlaps/caQTL_eQTL_pairs_betas_FC2.rds")
sign_count = dplyr::filter(priming_pairs, stimulation_state == "Stimulated") %>% 
  dplyr::select(gene_id, peak_id, snp_id, phenotype, beta, max_effect) %>% 
  dplyr::mutate(phenotype = ifelse(phenotype == "ATAC-seq", "ATAC", "RNA")) %>% 
  dplyr::distinct() %>% tidyr::spread(phenotype, beta) %>%
  dplyr::mutate(diff_sign = sign(ATAC)-sign(RNA))
