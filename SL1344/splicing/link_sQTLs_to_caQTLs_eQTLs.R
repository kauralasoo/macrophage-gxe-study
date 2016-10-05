library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")

#Import variant information
snp_info = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
snp_coords = dplyr::select(snp_info, snp_id, chr, pos)

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import eQTLs from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_qtl_df = extractQTLsFromList(rasqual_min_pvalues, fdr_cutoff = 0.1)
rasqual_qtl_coords = addVariantCoords(rasqual_qtl_df, vcf_file$snpspos)

#Import ATAC QTL variants
atac_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
atac_qtl_df = extractQTLsFromList(atac_min_pvalues, fdr_cutoff = 0.1)
atac_qtl_coords = addVariantCoords(atac_qtl_df, vcf_file$snpspos)

#Import Salmon sQTLs from disk
salmon_qtl_hits = readRDS("results/SL1344/salmon/salmon_qtl_hits.rds")
salmon_qtl_df = purrr::map_df(salmon_qtl_hits, identity, .id = "condition_name")

salmon_ensembl_qtl_hits = readRDS("results/SL1344/salmon/salmon_ensembl_qtl_hits.rds")
salmon_ensembl_qtl_df = purrr::map_df(salmon_ensembl_qtl_hits, identity, .id = "condition_name")


#Find overlaps
rna_atac_overlaps = findGWASOverlaps(filtered_pairs, atac_trait_pairs, vcf_file, max_distance = 5e5, min_r2 = 0.8)

sqtl_df = dplyr::filter(salmon_qtl_df, condition_name == "IFNg", event_type == "upstream") %>% 
  dplyr::select(gene_id, snp_id)
sqtl_df_contained = dplyr::filter(salmon_qtl_df, condition_name == "IFNg", event_type == "contained") %>% 
  dplyr::select(gene_id, snp_id)
sqtl_df_downstream = dplyr::filter(salmon_qtl_df, condition_name == "IFNg", event_type == "downstream") %>% 
  dplyr::select(gene_id, snp_id)

eqtl_df = dplyr::filter(rasqual_qtl_coords, condition_name == "IFNg") %>% dplyr::select(snp_id, chr, pos)
caqtl_df = dplyr::filter(atac_qtl_coords, condition_name == "IFNg") %>% dplyr::select(snp_id, chr, pos)

up_olap = (findGWASOverlaps(sqtl_df, eqtl_df, vcf_file, max_distance = 5e5, min_r2 = 0.9)$gene_id %>% 
  unique() %>% length())/885
up_olap_ca = (findGWASOverlaps(sqtl_df, caqtl_df, vcf_file, max_distance = 5e5, min_r2 = 0.9)$gene_id %>% 
                unique() %>% length())/885

contained_olap = (findGWASOverlaps(sqtl_df_contained, eqtl_df, vcf_file, max_distance = 5e5, min_r2 = 0.9)$gene_id %>% 
                    unique() %>% length())/560
contained_olap_ca = (findGWASOverlaps(sqtl_df_contained, caqtl_df, vcf_file, max_distance = 5e5, min_r2 = 0.9)$gene_id %>% 
                       unique() %>% length())/560

down_olap = (findGWASOverlaps(sqtl_df_downstream, eqtl_df, vcf_file, max_distance = 5e5, min_r2 = 0.9)$gene_id %>% 
               unique() %>% length())/871
down_olap_ca = (findGWASOverlaps(sqtl_df_downstream, caqtl_df, vcf_file, max_distance = 5e5, min_r2 = 0.9)$gene_id %>% 
                  unique() %>% length())/871

#Conclusion: it does not look very convincing








