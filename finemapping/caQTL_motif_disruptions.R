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
#Find indels
snp_info = dplyr::mutate(snp_info, is_indel = ifelse(pmax(nchar(alt), nchar(ref)) > 1, TRUE, FALSE))

#Keep only motfis that are enriched in macrophages
mf_enriched_motifs = read.table("results/ATAC/motif_analysis/cisBP_selected_enriched_motifs.txt", header = TRUE, stringsAsFactors = FALSE)
motif_names = dplyr::select(mf_enriched_motifs, motif_id, tf_name)
cisbp_pwm_enriched = cisbp_pwm_list[mf_enriched_motifs$motif_id]

#Import ATAC peak sequences from disk
sequences = Biostrings::readDNAStringSet("annotations/ATAC_consensus_peaks.fasta")
peak_ids = strsplit(names(sequences), "=") %>% lapply(function(x) x[2]) %>% unlist()
names(sequences) = peak_ids

#Extract IFNg-specific caQTLs
qtls = dplyr::filter(variable_qtls$appear, cluster_id == 3) %>% dplyr::select(gene_id, snp_id) %>% unique()

#Identify potential causal snps
cset_filtered = dplyr::filter(unique_peaks$lead_credible_sets, gene_id %in% qtls$gene_id)
overlap_genes = dplyr::semi_join(cset_filtered, qtls, by = c("gene_id","snp_id"))
credible_snps = dplyr::semi_join(unique_peaks$peak_snp_pairs, overlap_genes, by = "gene_id") %>%
  dplyr::filter(snp_count <= 2)

credible_snps_no_indel = dplyr::semi_join(credible_snps, dplyr::filter(snp_info, snp_id %in% credible_snps$snp_id) %>% 
                                            dplyr::filter(is_indel == FALSE), by = "snp_id")

#Calculate motif disruptions
motif_disruptions = quantifyMultipleVariants(credible_snps_no_indel, cisbp_pwm_enriched, atac_data$gene_metadata, sequences, snp_info)

#Filter the results
disruption_table_2 = dplyr::left_join(motif_disruptions, motif_names, by = "motif_id") %>%
  dplyr::group_by(gene_id, snp_id, motif_id) %>% 
  dplyr::arrange(-abs(rel_diff)) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-abs(rel_diff)) %>%
  dplyr::filter(rel_diff > 0.05)


#Indels
cs = dplyr::filter(credible_snps, gene_id == "ATAC_peak_256699")
motif_disruptions = quantifyMultipleVariants(cs, cisbp_pwm_enriched, atac_data$gene_metadata, sequences, snp_info)

d = quantifyMotifDisruption(cisbp_pwm_enriched[[1]], "ATAC_peak_256699", "rs75758825", atac_data$gene_metadata, sequences, snp_info)


a = quantifyMultipleMotifs("ATAC_peak_110001","rs11080327",cisbp_pwm_enriched, atac_data$gene_metadata, sequences, snp_info)

