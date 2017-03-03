library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")
library("GenomicRanges")

#Functions
testEnrichment <- function(fg_peaks, bg_peaks, fimo_hits, peak_metadata, unique_motifs, tf_name_casual){
  
  #Combine all peaks
  all_peaks = dplyr::bind_rows(fg_peaks, bg_peaks) %>% unique()
  
  #Calculate enrichment
  enrichment_result = fimoRelativeEnrichment(fg_peaks, all_peaks, fimo_hits, peak_metadata) %>%
    dplyr::left_join(unique_motifs, by = "motif_id") %>% 
    dplyr::arrange(fisher_pvalue) %>% 
    dplyr::filter(tf_name %in% c("RELA", "FOS", "IRF1", "SPI1")) %>% 
    dplyr::select(OR_log2, ci_lower_log2, ci_higher_log2, fisher_pvalue, tf_name) %>% 
    dplyr::left_join(tf_name_casual, by = "tf_name")
  return(enrichment_result)
}

#Load gene-peak pairs
pairs = readRDS("results/ATAC_RNA_overlaps/condition_specific_pairs.rds")
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Import all shared eQTLs-caQTLs
rna_atac_overlaps = readRDS("results/ATAC_RNA_overlaps/QTL_overlap_list_R2.rds")

#Import ATAC interaction test p-values
interaction_df = readRDS("results/ATAC/QTLs/rasqual_interaction_results.rds")
no_interaction_peaks = dplyr::filter(interaction_df, p_fdr > 0.1)

#Import RNA interaction test p-values
rna_interaction_df = readRDS("results/SL1344/eQTLs/SL1344_interaction_pvalues_lme4.rds")
no_interaction_genes = dplyr::filter(rna_interaction_df, p_fdr > 0.5)
no_interaction_pairs = dplyr::semi_join(rna_atac_overlaps, no_interaction_genes, by = "gene_id")

#Import condition-specific QTLs
var_qtls = readRDS("results/SL1344/eQTLs/appeat_disappear_eQTLs.rds")

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

#Are persistent peaks for condition-specific eQTLs enriched for condition-specific TF motifs

#Import FIMO motif matches
fimo_hits = readr::read_delim("results/ATAC/cisBP/FIMO_CISBP_results.long.txt.gz", delim = "\t", col_types = c("cciicddcc"), 
                              col_names = c("motif_id","seq_name","start","end","strand","score","p_value","dummy","matched_seq"), skip = 1)
fimo_hits_clean = tidyr::separate(fimo_hits, seq_name, c("prefix","gene_id"), sep = "=")

#Import motfis with lower p-value threshold
fimo_hits = readr::read_delim("results/ATAC/cisBP/FIMO_CISBP_results.1e-4.txt.gz", delim = "\t", col_types = c("cciicddcc"), 
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
  unique() 

#Add friendly names for motifs
interesting_tfs = c("RELA","IRF1","FOS","SPI1")
tf_name_casual = data_frame(new_name = c("NF-kB","IRF","AP-1","PU.1"), tf_name = interesting_tfs) %>%
  dplyr::mutate(new_name = factor(new_name, levels = rev(new_name)))

ifng_genes = unique(dplyr::filter(var_qtls$appear, new_cluster_id %in% c(5,6))$gene_id)
#Extract peaks that are present already in the naive state
ifng_peaks_present = dplyr::filter(pairs$IFNg, phenotype == "ATAC", condition_name == "naive", beta > 0.59)  %>% 
  dplyr::filter(gene_id %in% ifng_genes) %>%
  dplyr::select(peak_id) %>%
  dplyr::rename(gene_id = peak_id)
sl1344_peaks_present = dplyr::filter(pairs$SL1344, phenotype == "ATAC", condition_name == "naive", beta > 0.59)   %>% 
  dplyr::select(peak_id) %>%
  dplyr::rename(gene_id = peak_id)
ifng_sl1344_peaks_present = dplyr::filter(pairs$IFNg_SL1344, phenotype == "ATAC", condition_name == "naive", beta > 0.59)   %>% 
  dplyr::select(peak_id) %>%
  dplyr::rename(gene_id = peak_id)

#Calculate all of the enrichments
sl1344_enrichments = testEnrichment(sl1344_peaks_present, bg_peaks, fimo_hits_clean, atac_list$gene_metadata, unique_motifs, tf_name_casual) %>%
  dplyr::mutate(condition_name = "S")
ifng_enrichments = testEnrichment(ifng_peaks_present, bg_peaks, fimo_hits_clean, atac_list$gene_metadata, unique_motifs, tf_name_casual) %>%
  dplyr::mutate(condition_name = "I")
ifng_sl1344_enrichments = testEnrichment(ifng_sl1344_peaks_present, bg_peaks, fimo_hits_clean, atac_list$gene_metadata, unique_motifs, tf_name_casual) %>%
  dplyr::mutate(condition_name = "I+S")

#Make an enrichment plot
all_enrichment = dplyr::bind_rows(sl1344_enrichments,ifng_enrichments,ifng_sl1344_enrichments) %>% 
  dplyr::mutate(condition_name = factor(condition_name, levels = c("I","S","I+S")))

enrichment_plot = ggplot(all_enrichment, aes(y = new_name, x = OR_log2, xmin = ci_lower_log2, xmax = ci_higher_log2)) + 
  geom_point() + 
  geom_errorbarh(aes(height = 0)) +
  facet_wrap(~condition_name) + 
  xlab(expression(paste(Log[2], " fold enrichment", sep = ""))) +
  ylab("TF motif name") + 
  theme_light() +
  scale_x_continuous(expand = c(0, 0), limits = c(-3,3)) +
  theme(legend.key = element_blank()) + 
  geom_vline(aes(xintercept = 0), size = 0.3)
ggsave("figures/main_figures/caQTL_primed_motfis.pdf", plot = enrichment_plot, width = 5, height = 2.5)



#### Check for enrichment against ChIP-Seq peaks ####

#Import ChIP peaks from disk
naive_peaks = readRDS("results/ATAC/ChIP_enrichment/naive_combined_peaks.rds")
qiao_peaks = readRDS("results/ATAC/ChIP_enrichment/ivashkiv_peaks.rds")
wong_peaks = readRDS("results/ATAC/ChIP_enrichment/CIITA-RFX5_joint_peaks.rds")

#Extract peak coordinates
sl1344_persistent_ranges = dplyr::filter(atac_list$gene_metadata, gene_id %in% sl1344_persistent_peaks$gene_id) %>% 
  dplyr::transmute(seqnames = chr, start, end, strand = "+", gene_id) %>% dataFrameToGRanges()


persistent_ranges = dplyr::filter(atac_list$gene_metadata, gene_id %in% bg_peaks$gene_id) %>% 
  dplyr::transmute(seqnames = chr, start, end, strand = "+", gene_id) %>% dataFrameToGRanges()
ifng_ranges = dplyr::filter(atac_list$gene_metadata, gene_id %in% ifng_peaks_present$gene_id) %>% 
  dplyr::transmute(seqnames = chr, start, end, strand = "+", gene_id) %>% dataFrameToGRanges()
sl1344_ranges = dplyr::filter(atac_list$gene_metadata, gene_id %in% sl1344_peaks_present$gene_id) %>% 
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
