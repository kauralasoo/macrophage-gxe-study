library("Biostrings")
library("TFBSTools")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("purrr")

#Import peak counts
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")
atac_data$sample_metadata$condition_name = factor(atac_data$sample_metadata$condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344"))

#Import SNP coords and alleles
snp_info = readr::read_delim("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz", 
                               delim = "\t", col_types = "cdccc", col_names = c("chr","pos","snp_id","ref","alt"))

#Perform analysis with TFBSTools
cisbp_pwm_list = readRDS("results/ATAC/cisBP_PWMatrixList.rds")

#Import motif metadata
TF_information = readr::read_tsv("~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/TF_Information.txt")
colnames(TF_information)[6] = "gene_id"
unique_motifs = dplyr::select(TF_information, Motif_ID, gene_id, TF_Name) %>% dplyr::filter(Motif_ID != ".") %>%
  dplyr::rename(motif_id = Motif_ID, tf_name = TF_Name)

#Filter motifs by TF expression
combined_expression_data = readRDS("../macrophage-gxe-study/results/SL1344/combined_expression_data.rds")
mean_tpm = calculateMean(combined_expression_data$tpm, as.data.frame(combined_expression_data$sample_metadata), "condition_name")
expressed_genes = names(which(apply(mean_tpm, 1, max) > 1))
expressed_motifs = dplyr::filter(unique_motifs, gene_id %in% expressed_genes)

#Expressed_motifs
cisbp_pwm_expressed = cisbp_pwm_list[intersect(names(cisbp_pwm_list), expressed_motifs$motif_id)]

#Import ATAC peak sequences from disk
sequences = readDNAStringSet("annotations/ATAC_consensus_peaks.fasta")
peak_ids = strsplit(names(sequences), "=") %>% lapply(function(x) x[2]) %>% unlist()
names(sequences) = peak_ids

res = quantifyMotifDisruption(cisbp_pwm_list[["M6484_1.02"]], "ATAC_peak_145162", "rs7594476", atac_data$gene_metadata, sequences, snp_info) %>%
  dplyr::filter(abs(rel_diff) > 0, max_rel_score > 0.8)


#Find all possible disruptions
motif_disruptions = purrr::map(as.list(cisbp_pwm_expressed), ~quantifyMotifDisruption(., "ATAC_peak_145162", "rs7594476", atac_data$gene_metadata, sequences, snp_info) %>%
                                 dplyr::filter(abs(rel_diff) > 0, max_rel_score > 0.8))

dispruption_table = plyr::ldply(motif_disruptions)

filtered = dplyr::left_join(dispruption_table, expressed_motifs) %>% 
  dplyr::arrange(-abs(rel_diff)) %>%
  dplyr::group_by(tf_name) %>% 
  dplyr::filter(row_number() == 1)



