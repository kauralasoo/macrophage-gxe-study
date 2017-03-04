library("devtools")
library("plyr")
library("dplyr")
library("purrr")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#### Link eQTLs to caQTLs based on R2 overlap ####
#Import lead eQTL SNPs
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_qtl_df = purrr::map_df(rasqual_min_pvalues, 
                               ~dplyr::filter(., p_eigen < fdr_thresh), .id = "condition_name") %>% 
  dplyr::arrange(p_nominal)
joint_pairs = dplyr::select(rasqual_qtl_df, gene_id, snp_id) %>% unique() 
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

#Import ATAC QTL variants
atac_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
atac_qtl_df = purrr::map_df(atac_min_pvalues, 
                            ~dplyr::filter(., p_eigen < fdr_thresh), .id = "condition_name") %>% 
  dplyr::arrange(p_nominal)
atac_joint_pairs = dplyr::select(atac_qtl_df, gene_id, snp_id) %>% unique() 
atac_filtered_pairs = filterHitsR2(atac_joint_pairs, vcf_file$genotypes, .8)


#Find overlaps using the GWAS overlap code
atac_trait_pairs = addVariantCoords(atac_filtered_pairs, vcf_file$snpspos) %>%
  dplyr::rename(peak_id = gene_id)
rna_atac_overlaps = findGWASOverlaps(filtered_pairs, atac_trait_pairs, vcf_file, max_distance = 5e5, min_r2 = 0.8)
saveRDS(rna_atac_overlaps, "results/ATAC_RNA_overlaps/QTL_overlap_list_R2.rds")

#Identify shared QTLs
rna_atac_overlaps = readRDS("results/ATAC_RNA_overlaps/QTL_overlap_list_R2.rds")
shared_qtls = dplyr::select(rna_atac_overlaps, gene_id, peak_id) %>% unique()

#Find minimal p-values for each peaks across conditions
atac_unique_pvalues = purrr::map_df(atac_min_pvalues, identity, .id = "condition_name") %>%
  dplyr::filter(gene_id %in% unique(rna_atac_overlaps$peak_id)) %>%
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(peak_id = gene_id, p_nominal)

#Find unique pairs between genes and peaks
unique_pairs_r2 = dplyr::left_join(rna_atac_overlaps, atac_unique_pvalues, by = "peak_id") %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(gene_id, -R2, p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::filter(chr != "X") %>%
  dplyr::ungroup()

unique_pairs_minp = dplyr::left_join(rna_atac_overlaps, atac_unique_pvalues, by = "peak_id") %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(gene_id, p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::filter(chr != "X") %>%
  dplyr::ungroup()

#Save QTL pairs to disk
pair_list = list(R2 = unique_pairs_r2, minp = unique_pairs_minp)
saveRDS(pair_list, "results/ATAC_RNA_overlaps/ATAC_RNA_pairs.rds")




