library("devtools")
library("cqn")
library("dplyr")
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
  dplyr::rename(chr = seqnames)

#Construct peak metadata
peak_metadata = dplyr::select(atac_counts, gene_id, length) %>%
  dplyr::left_join(gc_content, by = "gene_id") %>% 
  dplyr::left_join(peak_coords, by = "gene_id") %>%
  dplyr::mutate(gene_name = gene_id)

#Extract count matrix
counts = dplyr::select(atac_counts, -gene_id, -length)
rownames(counts) = atac_counts$gene_id

#Use CQN to normalize the counts
atac_cqn = calculateCQN(counts, peak_metadata)

#Extract donor to genotype mapping
line_metadata = readRDS("../macrophage-gxe-study/macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
donor_geno_map = dplyr::select(line_metadata, donor, genotype_id) %>% unique()

#Construct sample metadata for atac
design_matrix = constructDesignMatrix_ATAC(colnames(counts))
atac_sample_meta = dplyr::left_join(design_matrix, donor_geno_map, by = "donor") %>% 
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))

#Construct a design matrix from the sample names

results_list = list(
  exprs_counts = counts,
  exprs_cqn = atac_cqn,
  design = design_matrix,
  gene_metadata = peak_metadata)

#Save processed data to disk
saveRDS(results_list, "results/ATAC/ATAC_combined_accessibility_data.rds")
