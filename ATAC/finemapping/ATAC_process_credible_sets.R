#Import condition-specific ATAC QTLs from disk
variable_qtls = readRDS("results/ATAC/QTLs/rasqual_appear_disappear_qtls.rds")

#Extract IFNg-specific caQTLs
qtls = dplyr::filter(variable_qtls$appear, cluster_id == 2) %>% dplyr::select(gene_id, snp_id) %>% unique()

#Fetch credible sets for each QTL
peak_cs = fetchCredibleSets(qtls, atac_list$gene_metadata, atac_list$gene_metadata, atac_tabix_list$IFNg, vcf_file, cis_window = 5e4)

#Calculate summary stats on the credible sets
cs_summary = summariseCredibleSets(peak_cs)

#Identify potential causal SNPs
potential_causal_snps = dplyr::filter(peak_cs, gene_id == overlap_peak_id)
