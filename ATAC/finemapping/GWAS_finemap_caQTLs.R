
#Import caQTLs
caqtl_coloc_hits = readRDS("results/SL1344/coloc/caQTL_coloc_posterior_hits.rds") %>%
  dplyr::select(gene_id, snp_id, trait, gwas_lead, gene_name) %>% unique() 


#Import unique caQTL peaks
unique_peaks = readRDS("results/ATAC/QTLs/unique_qtl_peaks.rds") 
clustered_peaks = readRDS("results/ATAC/QTLs/clustered_qtl_peaks.rds")

#Count overlaps
unique_olaps = dplyr::left_join(caqtl_coloc_hits, unique_peaks$lead_snps, by = "gene_id") %>% dplyr::filter(!is.na(snp_count))
cluster_olaps = dplyr::left_join(caqtl_coloc_hits, unique_peaks$lead_snps, by = "gene_id")