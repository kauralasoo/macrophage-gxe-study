#Scan all caQTLs for motif disruptions
library("devtools")
library("dplyr")
library("purrr")
load_all("../seqUtils/")
library("TFBSTools")

batch_id = NULL

####### Get batch id from disk ######
f <- file("stdin")
open(f)
batch_id = readLines(f) %>% as.numeric()
close(f)
####### END #######

#Import ATAC data
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")

#Import cisbp motifs
cisbp_pwm_list = readRDS("results/ATAC/cisBP_PWMatrixList.rds")

#Import SNP coords and alleles
snp_info = readr::read_delim("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz", 
                             delim = "\t", col_types = "cdccc", col_names = c("chr","pos","snp_id","ref","alt"))
#Find indels
snp_info = dplyr::mutate(snp_info, indel_length = pmax(nchar(alt), nchar(ref))) %>%
  dplyr::mutate(is_indel = ifelse(indel_length > 1, TRUE, FALSE))

#Keep only motfis that are enriched in macrophages
mf_enriched_motifs = read.table("results/ATAC/motif_analysis/cisBP_selected_enriched_motifs.txt", header = TRUE, stringsAsFactors = FALSE)
motif_names = dplyr::select(mf_enriched_motifs, motif_id, tf_name)
cisbp_pwm_enriched = cisbp_pwm_list[mf_enriched_motifs$motif_id]

#Import putative causal variants
unique_peaks = readRDS("results/ATAC/QTLs/unique_qtl_peaks.rds")
unique_peaks_filtered = dplyr::filter(unique_peaks$peak_snp_pairs, snp_count <= 3)

#### Split genes into batches ####
batches = splitIntoBatches(nrow(unique_peaks_filtered), 50)
if(!is.null(batch_id)){
  selection = batches == batch_id
  unique_peaks_filtered = unique_peaks_filtered[selection,]
}

#Import ATAC peak sequences from disk
sequences = Biostrings::readDNAStringSet("annotations/ATAC_consensus_peaks.fasta")
peak_ids = strsplit(names(sequences), "=") %>% lapply(function(x) x[2]) %>% unlist()
names(sequences) = peak_ids

#Calculate motif disruptions
motif_disruptions = quantifyMultipleVariants(unique_peaks_filtered, cisbp_pwm_enriched, atac_data$gene_metadata, sequences, snp_info)
#saveRDS(motif_disruptions, "results/ATAC/motif_analysis/unique_peaks_disrupted_motifs.rds")

#Save output from each batch
if(!is.null(batch_id)){
  output_file = file.path("results/ATAC/motif_analysis/", paste0("motif_disruption_batch_",batch_id, ".txt"))
  write.table(motif_disruptions, output_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}
                        

#Filter the results
#disruption_table_3 = dplyr::left_join(motif_disruptions, motif_names, by = "motif_id") %>%
#  dplyr::group_by(gene_id, snp_id, motif_id) %>% 
#  dplyr::arrange(-abs(rel_diff)) %>%
#  dplyr::filter(row_number() == 1) %>%
#  dplyr::ungroup() %>%
#  dplyr::arrange(-abs(rel_diff)) %>%
#  dplyr::filter(rel_diff > 0.05)


#Perform motif disruption analysis for a single peak and motif
res = quantifyMotifDisruption(cisbp_pwm_enriched[["M6119_1.02"]], "ATAC_peak_145162","rs7594476", atac_data$gene_metadata, sequences, snp_info)

