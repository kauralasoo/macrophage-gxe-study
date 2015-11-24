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
  dplyr::rename(chr = seqnames, left = start, right = end) %>%
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
atac_tpm = calculateTPM(counts, peak_metadata)
atac_tpm = atac_tpm[peak_metadata$gene_id,]

#Extract donor to genotype mapping
line_metadata = readRDS("../macrophage-gxe-study/macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
donor_geno_map = dplyr::select(line_metadata, donor, genotype_id) %>% unique()

#Construct sample metadata for atac
design_matrix = constructDesignMatrix_ATAC(colnames(counts))
atac_sample_meta = dplyr::left_join(design_matrix, donor_geno_map, by = "donor") %>% 
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))

#Construct separate design matrices
cond_A_design = dplyr::filter(atac_sample_meta, condition_char == "A")
cond_B_design = dplyr::filter(atac_sample_meta, condition_char == "B")
cond_C_design = dplyr::filter(atac_sample_meta, condition_char == "C")
cond_D_design = dplyr::filter(atac_sample_meta, condition_char == "D")
design_list = list(naive = cond_A_design, IFNg = cond_B_design, SL1344 = cond_C_design, IFNg_SL1344 = cond_D_design)

#Set up expression data for each condition
condA_exp = extractSubset(dplyr::filter(atac_sample_meta, condition_char == "A"), atac_cqn)
condB_exp = extractSubset(dplyr::filter(atac_sample_meta, condition_char == "B"), atac_cqn)
condC_exp = extractSubset(dplyr::filter(atac_sample_meta, condition_char == "C"), atac_cqn)
condD_exp = extractSubset(dplyr::filter(atac_sample_meta, condition_char == "D"), atac_cqn)
exprs_cqn_list = list(naive = condA_exp, IFNg = condB_exp, SL1344 = condC_exp, IFNg_SL1344 = condD_exp)

#Set up expression data for each condition
condA_tpm = extractSubset(dplyr::filter(atac_sample_meta, condition_char == "A"), atac_tpm)
condB_tpm = extractSubset(dplyr::filter(atac_sample_meta, condition_char == "B"), atac_tpm)
condC_tpm = extractSubset(dplyr::filter(atac_sample_meta, condition_char == "C"), atac_tpm)
condD_tpm = extractSubset(dplyr::filter(atac_sample_meta, condition_char == "D"), atac_tpm)
exprs_tpm_list = list(naive = condA_tpm, IFNg = condB_tpm, SL1344 = condC_tpm, IFNg_SL1344 = condD_tpm)

#Combine everything into a list
results_list = list(
  exprs_counts = counts,
  exprs_cqn = atac_cqn,
  exprs_tpm = atac_tpm,
  exprs_cqn_list = exprs_cqn_list,
  exprs_tpm_list = exprs_tpm_list,
  sample_meta = atac_sample_meta,
  design_list = design_list,
  gene_metadata = peak_metadata)

#Save processed data to disk
saveRDS(results_list, "results/ATAC/ATAC_combined_accessibility_data.rds")
