library("TFBSTools")

#Import QTL credible sets
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
credible_sets_df = importCredibleSets("results/ATAC/QTLs/rasqual_credible_sets.rds", atac_list$gene_metadata)

#Import caQTL and eQTL pairs
caQTL_eQTL_pairs = readRDS("results/ATAC_RNA_overlaps/caQTL_eQTL_pairs_betas.rds")
filtered_pairs = purrr::map_df(caQTL_eQTL_pairs, identity, .id = "condition") %>% 
  dplyr::filter(condition_name == "naive", phenotype == "ATAC-seq")

#Identify two groups of peaks
appear_peaks = dplyr::filter(filtered_pairs, abs(beta) < 0.59) %>% dplyr::select(peak_id) %>% unique()
persistent_peaks = dplyr::filter(filtered_pairs, abs(beta) > 0.59) %>% dplyr::select(peak_id) %>% unique()

#Identify credible sets of causal variants for both groups
appear_snps = purrr::map_df(credible_sets_df, ~dplyr::filter(., gene_id %in% appear_peaks$peak_id) %>% 
  dplyr::filter(gene_id == overlap_peak_id)) %>% 
  dplyr::select(gene_id, snp_id) %>% 
  unique() %>%
  dplyr::group_by(gene_id) %>% 
  dplyr::mutate(snp_count = length(snp_id)) %>%
  ungroup()

persistent_snps = purrr::map_df(credible_sets_df, ~dplyr::filter(., gene_id %in% persistent_peaks$peak_id) %>% 
                              dplyr::filter(gene_id == overlap_peak_id)) %>% 
  dplyr::select(gene_id, snp_id) %>% 
  unique() %>%
  dplyr::group_by(gene_id) %>% 
  dplyr::mutate(snp_count = length(snp_id)) %>%
  ungroup()

#Import cisbp motifs
cisbp_pwm_list = readRDS("results/ATAC/cisBP/cisBP_PWMatrixList.rds")

#Import SNP coords and alleles
snp_info = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

#Keep only motfis that are enriched in macrophages
mf_enriched_motifs = read.table("results/ATAC/motif_analysis/cisBP_selected_enriched_motifs.txt", header = TRUE, stringsAsFactors = FALSE)
motif_names = dplyr::select(mf_enriched_motifs, motif_id, tf_name)
cisbp_pwm_enriched = cisbp_pwm_list[mf_enriched_motifs$motif_id]

#Import ATAC peak sequences from disk
sequences = Biostrings::readDNAStringSet("annotations/chromatin/ATAC_consensus_peaks.fasta")
peak_ids = strsplit(names(sequences), "=") %>% lapply(function(x) x[2]) %>% unlist()
names(sequences) = peak_ids

#Calculate motif disruptions
persistent_disruptions = quantifyMultipleVariants(persistent_snps, cisbp_pwm_enriched, atac_list$gene_metadata, sequences, snp_info)
appear_disruptions = quantifyMultipleVariants(appear_snps, cisbp_pwm_enriched, atac_list$gene_metadata, sequences, snp_info)



