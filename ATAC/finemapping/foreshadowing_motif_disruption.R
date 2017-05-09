library("TFBSTools")
library("dplyr")

#Import QTL credible sets
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
credible_sets_df = importCredibleSets("results/ATAC/QTLs/rasqual_credible_sets.rds", atac_list$gene_metadata)

#Import caQTL and eQTL pairs
caQTL_eQTL_pairs = readRDS("results/ATAC_RNA_overlaps/caQTL_eQTL_pairs_betas.rds")
filtered_pairs = dplyr::filter(caQTL_eQTL_pairs, condition_name == "naive", phenotype == "ATAC-seq")

#Identify two groups of peaks
appear_peaks = dplyr::filter(filtered_pairs, abs(beta) < 0.59) %>% dplyr::select(peak_id) %>% unique()
persistent_peaks = dplyr::filter(filtered_pairs, abs(beta) > 0.59) %>% dplyr::select(peak_id) %>% unique()

#GBP1BA enhancer
persistent_peaks = data_frame(peak_id = "ATAC_peak_106417")

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

#Filter motif hits
persistent_hits = dplyr::left_join(persistent_disruptions, motif_names, by = "motif_id") %>% dplyr::filter(max_rel_score > 0.85, abs(rel_diff) > 0.03) %>% 
  dplyr::filter(snp_count <= 10) %>%
  dplyr::mutate(tf_group = ifelse(tf_name %in% c("FOS", "SPI1", "CEBPA", "CEBPB", "MAFB"), "Pioneer", ifelse(tf_name %in% c("IRF1","NFKB1","RELA"),"Stimulation", "Other")))

appear_hits = dplyr::left_join(appear_disruptions, motif_names, by = "motif_id") %>% dplyr::filter(max_rel_score > 0.85, abs(rel_diff) > 0.03) %>% 
  dplyr::filter(snp_count <= 10) %>%
  dplyr::mutate(tf_group = ifelse(tf_name %in% c("FOS", "SPI1", "CEBPA", "CEBPB", "MAFB"), "Pioneer", ifelse(tf_name %in% c("IRF1","NFKB1","RELA"),"Stimulation", "Other")))


#Count disruptions and 
persistent_f = dplyr::group_by(persistent_hits, tf_group, gene_id) %>% 
  dplyr::arrange(gene_id, tf_group, -rel_diff) %>% dplyr::filter(row_number() == 1) %>% ungroup()
appear_f = dplyr::group_by(appear_hits, tf_group, gene_id) %>% 
  dplyr::arrange(gene_id, tf_group, -rel_diff) %>% dplyr::filter(row_number() == 1) %>% ungroup()

table(persistent_f$tf_name)
length(unique(persistent_peaks$peak_id))

table(appear_f$tf_name)
length(unique(appear_peaks$peak_id))

fisher.test(matrix(c(9, 78-9, 0, 59), ncol = 2))


