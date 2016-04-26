library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Load p-values from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_min_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df = plyr::ldply(rasqual_min_hits, .id = "condition_name") %>% dplyr::arrange(p_nominal)
joint_pairs = dplyr::select(min_pvalue_df, gene_id, snp_id) %>% unique() 

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#List of RNA-seq tabix files
rna_list = list(naive = "results/SL1344/rasqual/output/naive_500kb/naive_500kb.sorted.txt.gz",
                IFNg = "results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.sorted.txt.gz",
                SL1344 = "results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.sorted.txt.gz",
                IFNg_SL1344 = "results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.sorted.txt.gz")
                

#List of ATAC tabix files
atac_tabix_list = list(naive = "../macrophage-chromatin/results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz",
                       IFNg = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz",
                       SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz",
                       IFNg_SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz")

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Fetch top genes
naive_list = idVectorToList(rasqual_min_hits$naive$gene_id)
naive_cs = purrr::map(naive_list, ~tabixFetchGenesQuick(.,rna_list$naive, combined_expression_data$gene_metadata, cis_window = 5e5)[[1]] %>%
                        dplyr::arrange(p_nominal) %>% addR2FromLead(vcf_file$genotypes) %>% dplyr::filter(R2 > 0.8))

#Find hits for naive QTLs
naive_annotated = purrr::map(naive_cs, ~annotateCredibleSet(.,atac_list$gene_metadata, atac_tabix_list$naive))
saveRDS(naive_annotated, "results/SL1344/eQTLs/naive_finemapped_QTLs.rds")
naive_filtered = purrr::map(naive_annotated, ~dplyr::filter(., overlap_peak_id == assoc_peak_id))

#Compare counts
raw_counts = lapply(naive_cs, nrow) %>% ldply(.id = "gene_id") %>% dplyr::rename(raw_counts = V1)
fm_counts = lapply(naive_filtered, nrow) %>% ldply(.id = "gene_id") %>% dplyr::rename(fm_counts = V1)
joint_counts = dplyr::left_join(raw_counts, fm_counts, by = "gene_id")


pairs = readRDS("results/SL1344/eQTLs/interaction_peak_gene_pairs.rds")


#Fetch top genes
ifng_list = idVectorToList(rasqual_min_hits$IFNg$gene_id)
ifng_list = ifng_list[intersect(pairs$IFNg$gene_id %>% unique(), rasqual_min_hits$IFNg$gene_id)]

#Construt credible sets for all eQTLs
ifng_cs = purrr::map(ifng_list, ~tabixFetchGenesQuick(.,rna_list$IFNg, combined_expression_data$gene_metadata, cis_window = 5e5)[[1]] %>%
                        dplyr::arrange(p_nominal) %>% dplyr::filter(chisq > 0.9*max(chisq)))
ifng_cs_df = purrr::map_df(ifng_cs, ~dplyr::mutate(.,chr = as.character(chr)))
ifng_cs_size = dplyr::group_by(ifng_cs_df, gene_id) %>% 
  summarise(snp_count = length(snp_id)) %>% 
  arrange(snp_count) %>% 
  dplyr::left_join(gene_name_map, by = "gene_id")

ifng_cs = purrr::map(ifng_list, ~tabixFetchGenesQuick(.,rna_list$IFNg, combined_expression_data$gene_metadata, cis_window = 5e5)[[1]] %>%
                       dplyr::arrange(p_nominal) %>%
                       dplyr::filter(cis_snp_r2 > 0.6) %>%
                       addR2FromLead(vcf_file$genotypes) %>% 
                       dplyr::filter(R2 > 0.8))
ifng_cs_df = purrr::map_df(ifng_cs, ~dplyr::mutate(.,chr = as.character(chr)))
ifng_cs_size2 = dplyr::group_by(ifng_cs_df, gene_id) %>% 
  summarise(snp_count = length(snp_id)) %>% 
  arrange(snp_count) %>% 
  dplyr::left_join(gene_name_map, by = "gene_id")


lead_snps = dplyr::group_by(ifng_cs_df, gene_id) %>% 
  dplyr::arrange(p_nominal) %>% dplyr::filter(row_number() == 1) %>%
  dplyr::left_join(ifng_cs_size2, by = "gene_id")

ifng_annotated = purrr::map(ifng_cs, ~annotateCredibleSet(.,atac_list$gene_metadata, atac_tabix_list$IFNg))
ifng_annotated_df = purrr::map_df(ifng_annotated, ~dplyr::mutate(.,chr = as.character(chr)))

olap_snp_count = dplyr::filter(ifng_annotated_df, !is.na(overlap_peak_id)) %>% 
  dplyr::group_by(gene_id) %>% dplyr::summarise(overlap_snp_count = length(snp_id))

caqtl_snp_count = dplyr::filter(ifng_annotated_df, overlap_peak_id == assoc_peak_id) %>% 
  dplyr::group_by(gene_id) %>% dplyr::summarise(caqtl_snp_count = length(snp_id))

olaps = dplyr::left_join(ifng_cs_size2, olap_snp_count) %>% dplyr::left_join(caqtl_snp_count)




