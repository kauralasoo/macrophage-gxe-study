library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")

#Load gene-peak pairs
pairs = readRDS("results/SL1344/eQTLs/interaction_peak_gene_pairs.rds")
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")

#IFNg
#Calculate diffs
ifng_atac_diff = dplyr::filter(pairs$IFNg, phenotype == "ATAC") %>% 
  group_by(gene_id) %>% 
  dplyr::arrange(condition_name) %>% 
  dplyr::mutate(scaled_diff = beta_scaled[2] - beta_scaled[1]) %>% 
  dplyr::mutate(diff = beta[2] - beta[1]) %>%
  dplyr::filter(condition_name == "naive") %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(scaled_diff)

sl1344_atac_diff = dplyr::filter(pairs$SL1344, phenotype == "ATAC") %>% 
  group_by(gene_id) %>% 
  dplyr::arrange(condition_name) %>% 
  dplyr::mutate(scaled_diff = beta_scaled[2] - beta_scaled[1]) %>% 
  dplyr::mutate(diff = beta[2] - beta[1]) %>%
  dplyr::filter(condition_name == "naive") %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(scaled_diff)

#Make lists of peaks
#IFNg
ifng_gained_peaks = dplyr::filter(ifng_atac_diff, diff > 0.32, scaled_diff > 0.2) %>% dplyr::select(peak_id) %>%
  dplyr::rename(gene_id = peak_id)
ifng_persistent_peaks = dplyr::filter(ifng_atac_diff, abs(diff) < 0.32) %>% dplyr::select(peak_id) %>%
  dplyr::rename(gene_id = peak_id)

#SL1344
sl1344_gained_peaks = dplyr::filter(sl1344_atac_diff, diff > 0.32, scaled_diff > 0.2) %>% dplyr::select(peak_id) %>%
  dplyr::rename(gene_id = peak_id)
sl1344_persistent_peaks = dplyr::filter(sl1344_atac_diff, abs(diff) < 0.32) %>% dplyr::select(peak_id) %>%
  dplyr::rename(gene_id = peak_id)

#Persistent
all_persistent = dplyr::bind_rows(sl1344_persistent_peaks, ifng_persistent_peaks) %>% unique()

#Merge all peaks together
all_peaks = dplyr::bind_rows(ifng_gained_peaks, ifng_persistent_peaks, sl1344_gained_peaks, sl1344_persistent_peaks) %>% unique()
all_all_peaks = dplyr::bind_rows(sl1344_atac_diff, ifng_atac_diff) %>% 
  dplyr::select(peak_id) %>% unique() %>% dplyr::rename(gene_id = peak_id)

#Import FIMO motif matches
fimo_hits = readr::read_delim("../macrophage-chromatin/results/ATAC/FIMO_CISBP_results.long.txt", delim = "\t", col_types = c("cciicddcc"), 
                              col_names = c("motif_id","seq_name","start","end","strand","score","p_value","dummy","matched_seq"), skip = 1)
fimo_hits_clean = tidyr::separate(fimo_hits, seq_name, c("prefix","gene_id"), sep = "=") 

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

#Enrichment
ifng_enrich = fimoRelativeEnrichment(ifng_gained_peaks, all_peaks, fimo_hits_clean, atac_list$gene_metadata)
ifng_gained_motifs = dplyr::left_join(ifng_enrich, unique_motifs) %>% 
  dplyr::semi_join(expressed_motifs, by = "motif_id") %>% 
  dplyr::mutate(p_fdr = p.adjust(p_hyper, method = "fdr")) %>%
  dplyr::arrange(p_hyper) %>% dplyr::filter(p_fdr < 0.2)
ifng_gained_motifs

sl1344_enrich = fimoRelativeEnrichment(sl1344_gained_peaks, all_peaks, fimo_hits_clean, atac_list$gene_metadata)
sl1344_gained_motifs = dplyr::left_join(sl1344_enrich, unique_motifs) %>% 
  dplyr::semi_join(expressed_motifs, by = "motif_id") %>% 
  dplyr::mutate(p_fdr = p.adjust(p_hyper, method = "fdr")) %>%
  dplyr::arrange(p_hyper) %>% dplyr::filter(p_fdr < 0.2)
sl1344_gained_motifs

persistent_enrich = fimoRelativeEnrichment(all_persistent, all_peaks, fimo_hits_clean, atac_list$gene_metadata)
persistent_gained_motifs = dplyr::left_join(persistent_enrich, unique_motifs) %>% 
  dplyr::semi_join(expressed_motifs, by = "motif_id") %>% 
  dplyr::mutate(p_fdr = p.adjust(p_hyper, method = "fdr")) %>%
  dplyr::arrange(p_hyper) %>% dplyr::filter(p_fdr < 0.2)
persistent_gained_motifs




#Identify motifs disrupted by variants

