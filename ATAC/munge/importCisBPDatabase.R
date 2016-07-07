library("dplyr")
library("DESeq2")
library("devtools")
library("edgeR")
library("rtracklayer")
load_all("macrophage-chromatin/housekeeping/")
load_all("../seqUtils/")
library("TFBSTools")

#Import TF information form disk
TF_information = readr::read_tsv("~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/TF_Information.txt")
colnames(TF_information)[6] = "gene_id"
unique_motifs = dplyr::select(TF_information, Motif_ID, gene_id, TF_Name) %>% dplyr::filter(Motif_ID != ".")

#Filter motifs by TF expression
#combined_expression_data = readRDS("../macrophage-gxe-study/results/SL1344/combined_expression_data.rds")
#mean_tpm = calculateMean(combined_expression_data$tpm, as.data.frame(combined_expression_data$sample_metadata), "condition_name")
#expressed_genes = names(which(apply(mean_tpm, 1, max) > 1))
#expressed_motifs = dplyr::filter(unique_motifs, gene_id %in% expressed_genes)

#Collapse motifs
motifs_collapsed = dplyr::group_by(unique_motifs, Motif_ID) %>% 
  dplyr::summarise(gene_id = paste(gene_id, collapse = ";"), TF_Name = paste(TF_Name,collapse = ";"), TF_count = length(Motif_ID))

#Import PWMs
cisbp_dir = "~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/pwms_all_motifs/"
pwmatrix_list = cisBPImportRecords(motifs_collapsed, cisbp_dir)
pfmatrix_list = cisBPImportRecords(motifs_collapsed, cisbp_dir, toPWM = FALSE)

#Export PWM list to disk
saveRDS(pwmatrix_list, "results/ATAC/cisBP/cisBP_PWMatrixList.rds")
saveRDS(pfmatrix_list, "results/ATAC/cisBP/cisBP_PFMatrixList.rds")
saveRDS(motifs_collapsed,"results/ATAC/cisBP/cisBP_motif_metadata.rds")

