load_all("../seqUtils/")

#Import ATAC data
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")

#Import condition-specific ATAC QTLs from disk
variable_qtls = readRDS("results/ATAC/QTLs/rasqual_appear_disappear_qtls.rds")

#Import putative causal variants
unique_peaks = readRDS("results/ATAC/QTLs/unique_qtl_peaks.rds")

#Import cisbp motifs
cisbp_pwm_list = readRDS("results/ATAC/cisBP_PWMatrixList.rds")

#Import SNP coords and alleles
snp_info = readr::read_delim("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz", 
                             delim = "\t", col_types = "cdccc", col_names = c("chr","pos","snp_id","ref","alt"))

#Keep only motfis that are enriched in macrophages
mf_enriched_motifs = read.table("results/ATAC/motif_analysis/cisBP_selected_enriched_motifs.txt", header = TRUE)
motif_names = dplyr::select(mf_enriched_motifs, motif_id, tf_name)
cisbp_pwm_enriched = cisbp_pwm_list[mf_enriched_motifs$motif_id]

#Import ATAC peak sequences from disk
sequences = Biostrings::readDNAStringSet("annotations/ATAC_consensus_peaks.fasta")
peak_ids = strsplit(names(sequences), "=") %>% lapply(function(x) x[2]) %>% unlist()
names(sequences) = peak_ids

#Extract IFNg-specific caQTLs
qtls = dplyr::filter(variable_qtls$appear, cluster_id == 2) %>% dplyr::select(gene_id, snp_id) %>% unique()

#Identify potential causal snps
cset_filtered = dplyr::filter(unique_peaks$lead_credible_sets, gene_id %in% qtls$gene_id)
overlap_genes = dplyr::semi_join(cset_filtered, qtls, by = c("gene_id","snp_id"))
credible_snps = dplyr::semi_join(unique_peaks$peak_snp_pairs, overlap_genes, by = "gene_id") %>%
  dplyr::filter(snp_count <= 2)

#Calculate motif disruptions
motif_disruptions = quantifyMultipleVariants(credible_snps, cisbp_pwm_enriched, atac_data$gene_metadata, sequences, snp_info)

#Filter the results
disruption_table = purrr::map_df(motif_disruptions, identity) %>% 
  dplyr::left_join(motif_names, by = "motif_id") %>%
  dplyr::arrange(-abs(rel_diff)) %>%
  dplyr::group_by(motif_id) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
