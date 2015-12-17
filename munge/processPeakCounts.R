library("devtools")
library("cqn")
library("dplyr")
library("rtracklayer")
load_all("../seqUtils/")
load_all("macrophage-chromatin/housekeeping/")

#Import raw peak counts
atac_counts = readRDS("results/ATAC/ATAC_combined_counts.rds")
gc_content = read.table("annotations/ATAC_Seq_joint_peaks.GC_content.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::rename(gene_id = X9_usercol, percentage_gc_content = X11_pct_gc)

#Import peak coordinates
peak_coords = rtracklayer::import.gff3("annotations/ATAC_Seq_joint_peaks.gff3") %>% 
  GenomicRanges::as.data.frame() %>% tbl_df() %>% 
  dplyr::select(gene_id, seqnames, start, end) %>% 
  dplyr::rename(chr = seqnames) %>%
  dplyr::mutate(strand = "+", score = 1000) 
  
#Construct peak metadata
peak_metadata = dplyr::select(atac_counts, gene_id, length) %>%
  dplyr::left_join(gc_content, by = "gene_id") %>% 
  dplyr::left_join(peak_coords, by = "gene_id") %>%
  dplyr::mutate(gene_name = gene_id) %>%
  dplyr::filter( !(chr %in% c("MT","Y")) )

#Extract count matrix
counts = dplyr::select(atac_counts, -gene_id, -length)
rownames(counts) = atac_counts$gene_id
counts = counts[peak_metadata$gene_id,]

#Use CQN to normalize the counts
atac_cqn = calculateCQN(counts, peak_metadata)
atac_cqn = atac_cqn[peak_metadata$gene_id,]

#Normalize data using TPM
atac_tpm = calculateTPM(counts, peak_metadata, fragment_length = 50)
atac_tpm = atac_tpm[peak_metadata$gene_id,]

#Extract donor to genotype mapping
atac_metadata = readRDS("macrophage-chromatin/data/SL1344/compiled_atac_metadata.rds")

#Combine everything into a list
results_list = list(
  counts = counts,
  cqn = atac_cqn,
  tpm = atac_tpm,
  sample_metadata = atac_metadata,
  gene_metadata = peak_metadata)

#Save processed data to disk
saveRDS(results_list, "results/ATAC/ATAC_combined_accessibility_data.rds")
