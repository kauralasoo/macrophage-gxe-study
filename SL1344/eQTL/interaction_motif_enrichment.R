library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")
library("GenomicRanges")

#Load gene-peak pairs
pairs = readRDS("results/ATAC_RNA_overlaps/condition_specific_pairs.rds")
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Import all shared eQTLs-caQTLs
rna_atac_overlaps = readRDS("results/ATAC_RNA_overlaps/QTL_overlap_list_R2.rds")

#Import ATAC interaction test p-values
interaction_df = readRDS("results/ATAC/QTLs/rasqual_interaction_results.rds")
no_interaction_peaks = dplyr::filter(interaction_df, p_fdr > 0.1)

#Import RNA interaction test p-values
rna_interaction_df = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues.rds")
no_interaction_genes = dplyr::filter(rna_interaction_df, p_fdr > 0.5)
no_interaction_pairs = dplyr::semi_join(rna_atac_overlaps, no_interaction_genes, by = "gene_id")

#Find minimal p-values
atac_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
peak_min_pvalues = purrr::map_df(atac_min_pvalues, ~dplyr::semi_join(., rna_atac_overlaps, by = c("gene_id" = "peak_id"))) %>% 
  dplyr::select(gene_id, p_nominal) %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::rename(peak_id = gene_id)

#Find the most associated peak for each gene
rna_atac_min_pairs = dplyr::left_join(rna_atac_overlaps, peak_min_pvalues, by = "peak_id") %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()

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
ifng_persistent_peaks = dplyr::filter(ifng_atac_diff, diff < 0.32) %>% dplyr::select(peak_id) %>%
  dplyr::rename(gene_id = peak_id)

#SL1344
sl1344_gained_peaks = dplyr::filter(sl1344_atac_diff, diff > 0.32, scaled_diff > 0.2) %>% dplyr::select(peak_id) %>%
  dplyr::rename(gene_id = peak_id)
sl1344_persistent_peaks = dplyr::filter(sl1344_atac_diff, diff < 0.32) %>% dplyr::select(peak_id) %>%
  dplyr::rename(gene_id = peak_id)

#Are persistent peaks for condition-specific eQTLs enriched for condition-specific TF motifs

#Import FIMO motif matches
fimo_hits = readr::read_delim("results/ATAC/cisBP/FIMO_CISBP_results.long.txt.gz", delim = "\t", col_types = c("cciicddcc"), 
                              col_names = c("motif_id","seq_name","start","end","strand","score","p_value","dummy","matched_seq"), skip = 1)
fimo_hits_clean = tidyr::separate(fimo_hits, seq_name, c("prefix","gene_id"), sep = "=") 

#Import motif metadata
TF_information = readr::read_tsv("~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/TF_Information.txt")
colnames(TF_information)[6] = "gene_id"
unique_motifs = dplyr::select(TF_information, Motif_ID, gene_id, TF_Name) %>% dplyr::filter(Motif_ID != ".") %>%
  dplyr::rename(motif_id = Motif_ID, tf_name = TF_Name)

#Calculate enrichments
bg_peaks = dplyr::select(rna_atac_min_pairs, peak_id) %>% 
  dplyr::semi_join(no_interaction_pairs, by = "peak_id") %>% 
  dplyr::rename(gene_id = peak_id) %>% 
  #dplyr::semi_join(no_interaction_peaks, by = "gene_id") %>% 
  unique() 
bg_peaks = dplyr::bind_rows(sl1344_persistent_peaks, ifng_persistent_peaks, bg_peaks) %>% unique()

sl1344_enrichment = fimoRelativeEnrichment(sl1344_persistent_peaks, bg_peaks, fimo_hits_clean, 
                                           atac_list$gene_metadata)

dplyr::left_join(sl1344_enrichment, unique_motifs, by = "motif_id") %>% 
  dplyr::arrange(fisher_pvalue) %>% 
  dplyr::filter(tf_name %in% c("RELA", "FOS", "IRF1", "SPI1", "STAT1")) %>% 
  dplyr::select(OR, ci_lower, ci_higher, fisher_pvalue, tf_name)

ifng_enrichment = fimoRelativeEnrichment(ifng_persistent_peaks, bg_peaks, fimo_hits_clean, 
                                         atac_list$gene_metadata)
dplyr::left_join(ifng_enrichment, unique_motifs, by = "motif_id") %>% 
  dplyr::arrange(fisher_pvalue) %>% 
  dplyr::filter(tf_name %in% c("RELA", "FOS", "IRF1", "SPI1", "STAT1")) %>% 
  dplyr::select(OR, ci_lower, ci_higher, fisher_pvalue, tf_name)



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

#Import list of enriched motifs
enriched_motifs = read.table("../macrophage-chromatin/results/ATAC/motif_analysis/cisBP_selected_enriched_motifs.txt", header = TRUE,
                             stringsAsFactors = FALSE)

#Import motif metadata
TF_information = readr::read_tsv("~/annotations/CisBP/Homo_sapiens_2016_03_10_11-59_am/TF_Information.txt")
colnames(TF_information)[6] = "gene_id"
unique_motifs = dplyr::select(TF_information, Motif_ID, gene_id, TF_Name) %>% dplyr::filter(Motif_ID != ".") %>%
  dplyr::rename(motif_id = Motif_ID, tf_name = TF_Name)

#Enrichment
ifng_enrich = fimoRelativeEnrichment(ifng_gained_peaks, all_peaks, fimo_hits_clean, atac_list$gene_metadata)
ifng_gained_motifs = dplyr::left_join(ifng_enrich, unique_motifs) %>% 
  dplyr::semi_join(enriched_motifs, by = "motif_id") %>% 
  dplyr::mutate(p_fdr = p.adjust(p_hyper, method = "fdr")) %>%
  dplyr::arrange(p_hyper) %>% dplyr::filter(p_fdr < 0.2)
ifng_gained_motifs

sl1344_enrich = fimoRelativeEnrichment(sl1344_gained_peaks, all_peaks, fimo_hits_clean, atac_list$gene_metadata)
sl1344_gained_motifs = dplyr::left_join(sl1344_enrich, unique_motifs) %>% 
  dplyr::semi_join(enriched_motifs, by = "motif_id") %>% 
  dplyr::mutate(p_fdr = p.adjust(p_hyper, method = "fdr")) %>%
  dplyr::arrange(p_hyper) %>% dplyr::filter(p_fdr < 0.2)
sl1344_gained_motifs

persistent_enrich = fimoRelativeEnrichment(all_persistent, all_peaks, fimo_hits_clean, atac_list$gene_metadata)
persistent_gained_motifs = dplyr::left_join(persistent_enrich, unique_motifs) %>% 
  dplyr::semi_join(unique_motifs, by = "motif_id") %>% 
  dplyr::mutate(p_fdr = p.adjust(p_hyper, method = "fdr")) %>%
  dplyr::arrange(p_hyper) %>% dplyr::filter(p_fdr < 0.2)
persistent_gained_motifs

#Identify motifs disrupted by variants


#### Check for enrichment against ChIP-Seq peaks ####

#Import ChIP peaks from disk
naive_peaks = readRDS("../macrophage-chromatin/results/ATAC/ChIP_enrichment/naive_combined_peaks.rds")
qiao_peaks = readRDS("../macrophage-chromatin/results/ATAC/ChIP_enrichment/ivashkiv_peaks.rds")
wong_peaks = readRDS("../macrophage-chromatin/results/ATAC/ChIP_enrichment/CIITA-RFX5_joint_peaks.rds")

#Extract peak coordinates
persistent_ranges = dplyr::filter(atac_list$gene_metadata, gene_id %in% all_persistent$gene_id) %>% 
  dplyr::transmute(seqnames = chr, start, end, strand = "+", gene_id) %>% dataFrameToGRanges()
sl1344_persistent_ranges = dplyr::filter(atac_list$gene_metadata, gene_id %in% sl1344_persistent_peaks$gene_id) %>% 
  dplyr::transmute(seqnames = chr, start, end, strand = "+", gene_id) %>% dataFrameToGRanges()
ifng_ranges = dplyr::filter(atac_list$gene_metadata, gene_id %in% ifng_gained_peaks$gene_id) %>% 
  dplyr::transmute(seqnames = chr, start, end, strand = "+", gene_id) %>% dataFrameToGRanges()
sl1344_ranges = dplyr::filter(atac_list$gene_metadata, gene_id %in% sl1344_gained_peaks$gene_id) %>% 
  dplyr::transmute(seqnames = chr, start, end, strand = "+", gene_id) %>% dataFrameToGRanges()
ranges_list = list(persistent = persistent_ranges, IFNg = ifng_ranges, SL1344 = sl1344_ranges)

#Count overlaps with peaks
stat1_olaps = purrr::map(ranges_list, ~GenomicRanges::findOverlaps(., qiao_peaks$STAT1_rep2_B) %>% 
                           subjectHits() %>% length()) %>% unlist()  %>% t() %>% as.data.frame()
irf1_olaps = purrr::map(ranges_list, ~GenomicRanges::findOverlaps(., qiao_peaks$IRF1_B) %>% 
                           subjectHits() %>% length()) %>% unlist()  %>% t() %>% as.data.frame()
pu1_olaps = purrr::map(ranges_list, ~GenomicRanges::findOverlaps(., naive_peaks[naive_peaks$name == "PU.1_Schmidt"]) %>% 
                          subjectHits() %>% length()) %>% unlist()  %>% t() %>% as.data.frame()
pu1_olaps_2 = purrr::map(ranges_list, ~GenomicRanges::findOverlaps(., naive_peaks[naive_peaks$name == "PU.1_Pham"]) %>% 
                         subjectHits() %>% length()) %>% unlist()  %>% t() %>% as.data.frame()
cebpb_olaps = purrr::map(ranges_list, ~GenomicRanges::findOverlaps(., naive_peaks[naive_peaks$name == "CEBPb_Reschen"]) %>% 
                           subjectHits() %>% length()) %>% unlist()  %>% t() %>% as.data.frame()
cebpb_olaps_2 = purrr::map(ranges_list, ~GenomicRanges::findOverlaps(., naive_peaks[naive_peaks$name == "CEBPb_Pham"]) %>% 
                           subjectHits() %>% length()) %>% unlist()  %>% t() %>% as.data.frame()

#Calculate enrichments of overlap
results = dplyr::bind_rows(stat1_olaps, irf1_olaps, pu1_olaps, pu1_olaps_2, cebpb_olaps, cebpb_olaps_2) %>% 
  dplyr::mutate(chip = c("STAT1_IFNg","IRF1_IFNg", "PU.1_Schmidt","PU.1_Pham","CEBPb_Reschen","CEBPb_Pham")) %>%
  dplyr::mutate(persistent_count = length(ranges_list$persistent), IFNg_count = length(ranges_list$IFNg), 
                SL1344_count = length(ranges_list$SL1344)) %>%
  dplyr::mutate(persistent_frac = persistent/persistent_count, IFNg_frac = IFNg/IFNg_count, SL1344_frac = SL1344/SL1344_count) %>%
  dplyr::group_by(chip) %>% 
  dplyr::mutate(IFNg_p = fisher.test(matrix(c(persistent_count-persistent,IFNg_count-IFNg,persistent,IFNg), nrow = 2))$p.value) %>%
  dplyr::mutate(SL1344_p = fisher.test(matrix(c(persistent_count-persistent,SL1344_count-SL1344,persistent,SL1344), nrow = 2))$p.value) %>%
  dplyr::select(chip, everything())

write.table(results, "results/SL1344/eQTLs/properties/condition_specific_TF_enrichments.txt", sep = "\t", quote = FALSE, row.names = FALSE)
