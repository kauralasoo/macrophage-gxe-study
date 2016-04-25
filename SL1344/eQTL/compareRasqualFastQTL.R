library("purrr")
library("devtools")
library("dplyr")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")

#Import p-value from both methods
rasqual_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
fastqtl_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")

#Identify genes tested by rasqual
gene_id_list = purrr::map(rasqual_pvalues, ~dplyr::select(., gene_id) %>% unlist())
intersect_genes = purrr::reduce(gene_id_list, intersect)
union_genes = purrr::reduce(gene_id_list, union)

#Keep only those fastQTl results that were also present in rasqual output
fastqtl_pvalues_filtered = purrr::map2(fastqtl_pvalues, rasqual_pvalues, function(x,y){
  test_df = dplyr::select(y, gene_id, n_tests)
  res = dplyr::filter(x, gene_id %in% y$gene_id) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::left_join(test_df, by = "gene_id") %>%
    dplyr::mutate(p_eigen = p.adjust(p_nominal, method = "bonferroni", n = n_tests)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_eigen_fdr = p.adjust(p_eigen, method = "fdr")) %>%
    dplyr::arrange(p_nominal)
  return(res)
  })

#Count the number of eQTLs detected with each method
fastqtl_qtl_count = map(fastqtl_pvalues_filtered, ~dplyr::filter(., p_eigen_fdr < 0.1) %>% nrow()) %>% unlist()
rasqual_qtl_count = map(rasqual_pvalues, ~dplyr::filter(., p_fdr < 0.1) %>% nrow()) %>% unlist()
qtl_counts = rbind(fastqtl_qtl_count, rasqual_qtl_count) %>% t() %>% as.data.frame()
qtl_counts = dplyr::mutate(qtl_counts, condition = rownames(qtl_counts)) %>% 
  dplyr::select(condition, everything()) %>% 
  dplyr::mutate(extra_qtls = (rasqual_qtl_count - fastqtl_qtl_count)/fastqtl_qtl_count)
write.table(qtl_counts, "results/SL1344/eQTLs/properties/fastQTL_vs_rasqual_eQTL_counts.txt", quote = FALSE, row.names = FALSE, sep ="\t")

#Join rasqual and fastqtl data frames for scatter plots
joint_pvalues = purrr::map2(rasqual_pvalues, fastqtl_pvalues_filtered, function(rasqual, fastqtl){
  p_rasqual = dplyr::transmute(rasqual, gene_id, p_rasqual = ifelse(p_nominal == 0, 1e-300, p_nominal)) %>%
    dplyr::mutate(p_rasqual_log = -log(p_rasqual, 10))
  p_fastqtl = dplyr::transmute(fastqtl, gene_id, gene_name, p_fastqtl = p_nominal) %>%
    dplyr::mutate(p_fastqtl_log = -log(p_fastqtl, 10))
  result = dplyr::left_join(p_rasqual, p_fastqtl, by = "gene_id")
})
joint_pvalues_df = map_df(joint_pvalues, identity, .id = "condition_name") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))

#Make scatter plots
label_genes = c("CHURC1", "PTK2B", "ERAP2")
label_genes = snp_overlaps$gene_name
scatter_plots = ggplot(joint_pvalues_df, aes(x = p_fastqtl_log, y = p_rasqual_log, label = gene_name)) + 
  geom_point() + 
  geom_text(data = dplyr::filter(joint_pvalues_df, gene_name %in% label_genes), color = "blue") +
  facet_wrap(~condition_name) +
  ylab("RASQUAL p-value (-log10)") +
  xlab("FastQTL p-value (-log10)")
ggsave("results/SL1344/eQTLs/properties/fastQTL_vs_rasqual_eQTL_scatters.pdf", plot = scatter_plots, width = 10, height = 10)



#### Identify genes that could be affected by genotyping errors ####
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import SNP coordinates
snp_info = readr::read_delim("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz", 
                             delim = "\t", col_types = "cdccc", col_names = c("chr","pos","snp_id","ref","alt"))
bad_snps = read.table("genotypes/SL1344/imputed_20151005/error_snps_bf20.txt")
bad_snp_info = dplyr::filter(snp_info, snp_id %in% bad_snps$V3)

snp_overlaps = countSnpsOverlapingExons(combined_expression_data$gene_metadata, bad_snp_info, cis_window = 500000) %>%
  dplyr::select(gene_id, feature_snp_count) %>% 
  dplyr::filter(feature_snp_count > 0) %>% 
  dplyr::left_join(gene_name_map, by = "gene_id")
write.table(snp_overlaps$gene_id, "macrophage-gxe-study/data/rasqual_bad_fSNP_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Identify bad feature SNPs
bad_fSNPs = countSnpsOverlapingExons(combined_expression_data$gene_metadata, bad_snp_info, return_fSNPs = TRUE)
bad_fSNP_ids = dplyr::filter(bad_snp_info, pos %in% start(bad_fSNPs))$snp_id
write.table(bad_fSNP_ids, "macrophage-gxe-study/data/rasqual_bad_feature_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)



