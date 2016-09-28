library("readr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("purrr")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")

#Import p-value from both methods
rasqual_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
fastqtl_pvalues = readRDS("results/ATAC/QTLs/fastqtl_min_pvalues.rds")

#Identify genes tested by rasqual
gene_id_list = purrr::map(rasqual_pvalues, ~dplyr::select(., gene_id) %>% unlist())
intersect_genes = purrr::reduce(gene_id_list, intersect)
union_genes = purrr::reduce(gene_id_list, union)

#Keep only those fastQTL results that were also present in rasqual output
fastqtl_pvalues_filtered = purrr::map2(fastqtl_pvalues, rasqual_pvalues, function(x,y){
  res = dplyr::filter(x, gene_id %in% y$gene_id) %>%
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
write.table(qtl_counts, "results/ATAC/QTLs/properties/fastQTL_vs_rasqual_QTL_counts.txt", quote = FALSE, row.names = FALSE, sep ="\t")

#Join rasqual and fastqtl data frames for scatter plots
joint_pvalues = purrr::map2(rasqual_pvalues, fastqtl_pvalues_filtered, function(rasqual, fastqtl){
  p_rasqual = dplyr::transmute(rasqual, gene_id, p_rasqual = ifelse(p_nominal == 0, 1e-300, p_nominal)) %>%
    dplyr::mutate(p_rasqual_log = -log(p_rasqual, 10))
  p_fastqtl = dplyr::transmute(fastqtl, gene_id, p_fastqtl = p_nominal) %>%
    dplyr::mutate(p_fastqtl_log = -log(p_fastqtl, 10))
  result = dplyr::left_join(p_rasqual, p_fastqtl, by = "gene_id")
})
joint_pvalues_df = map_df(joint_pvalues, identity, .id = "condition_name") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))


#Make scatter plots
scatter_plots = ggplot(joint_pvalues_df, aes(x = p_fastqtl_log, y = p_rasqual_log)) + 
  geom_point() + 
  facet_wrap(~condition_name) +
  ylab("RASQUAL p-value (-log10)") +
  xlab("FastQTL p-value (-log10)")
ggsave("results/ATAC/QTLs/properties/fastQTL_vs_rasqual_eQTL_scatters.pdf", plot = scatter_plots, width = 10, height = 10)


#Make Q-Q plots
#Calculate expected p-values
rasqual_expected = purrr::map(rasqual_pvalues, ~addExpectedPvalue(.))
fastqtl_expected = purrr::map(fastqtl_pvalues_filtered, ~addExpectedPvalue(.))

#Collect p-values for the Q-Q plots
qqplot_lists = purrr::map2(rasqual_expected, fastqtl_expected, function(rasqual, fastqtl){
  rasqual_res = dplyr::transmute(rasqual, p_eigen, p_expected, method = "RASQUAL")
  fastqtl_res = dplyr::transmute(fastqtl, p_eigen, p_expected, method = "FastQTL")
  res = rbind(rasqual_res, fastqtl_res)
})
qqplot_df = purrr::map_df(qqplot_lists, identity, .id = "condition_name") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344", "IFNg_SL1344")))

#Make a Q-Q plot for each condition
qqplot = ggplot(qqplot_df, aes(x = -log(p_expected,10), y = -log(p_eigen, 10),color = method)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  facet_wrap(~condition_name, ncol = 4) +
  theme_light() +
  ylab(expression(paste(Log[10], " observed p-value", sep = ""))) +
  xlab(expression(paste(Log[10], " expected p-value", sep = "")))
ggsave("figures/supplementary/fastQTL_vs_rasqual_caQTL_qqplots.png", plot = qqplot, width = 8, height = 4)
ggsave("results/ATAC/QTLs/properties/fastQTL_vs_rasqual_eQTL_qqplots.pdf", plot = qqplot, width = 10, height = 10)



#See what is the effect if we only look at peaks that contain at least one feature SNP
snp_info = readr::read_delim("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz", 
                             delim = "\t", col_types = "cdccc", col_names = c("chr","pos","snp_id","ref","alt"))
snp_info = dplyr::mutate(snp_info, indel_length = pmax(nchar(alt), nchar(ref))) %>%
  dplyr::mutate(is_indel = ifelse(indel_length > 1, TRUE, FALSE))

#Import data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
gene_name_map = dplyr::select(atac_list$gene_metadata, gene_id, gene_name)

#Test for overlap
peak_ranges = dplyr::transmute(atac_list$gene_metadata, gene_id = gene_id, seqnames = chr, start, end, strand = "+") %>% dataFrameToGRanges()
snp_ranges = dplyr::filter(snp_info, is_indel == FALSE) %>%
  dplyr::transmute(snp_id, seqnames = chr, start = pos, end = pos, strand = "+") %>% 
  dataFrameToGRanges()
peaks_with_snps = peak_ranges[GenomicRanges::queryHits(GenomicRanges::findOverlaps(peak_ranges, snp_ranges)) %>% unique()]$gene_id

#Keep only peaks with SNPs
rasqual_pvalues_snps = purrr::map(rasqual_pvalues, ~dplyr::filter(.,gene_id %in% peaks_with_snps))
fastqtl_pvalues_snps = purrr::map(fastqtl_pvalues_filtered, ~dplyr::filter(.,gene_id %in% peaks_with_snps))

#Count qtls
fastqtl_qtl_count = map(fastqtl_pvalues_snps, ~dplyr::filter(., p_eigen_fdr < 0.1) %>% nrow()) %>% unlist()
rasqual_qtl_count = map(rasqual_pvalues_snps, ~dplyr::filter(., p_fdr < 0.1) %>% nrow()) %>% unlist()
qtl_counts = rbind(fastqtl_qtl_count, rasqual_qtl_count) %>% t() %>% as.data.frame()
qtl_counts = dplyr::mutate(qtl_counts, condition = rownames(qtl_counts)) %>% 
  dplyr::select(condition, everything()) %>% 
  dplyr::mutate(extra_qtls = (rasqual_qtl_count - fastqtl_qtl_count)/fastqtl_qtl_count)
write.table(qtl_counts, "results/ATAC/QTLs/properties/fastQTL_vs_rasqual_QTL_counts_contain_snps.txt", 
            quote = FALSE, row.names = FALSE, sep ="\t")
