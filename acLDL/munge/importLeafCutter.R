library("dplyr")
library("SummarizedExperiment")
library("devtools")
load_all("../seqUtils/")

#Import Ensembl quant results for reference
ensembl_quants = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
sample_names = colnames(ensembl_quants)

#Import intron and cluster counts from leafcutter
intron_counts = read.table("processed/acLDL/leafcutter/leafcutter_perind_numers.counts.gz")[,sample_names]

#Extract intron coords from counts
intron_metadata = dplyr::data_frame(transcript_id = rownames(intron_counts)) %>%
  tidyr::separate(transcript_id, c("chr","transcript_start","transcript_end","gene_id"), sep = ":", remove = FALSE) %>%
  dplyr::mutate(transcript_start = as.integer(transcript_start), transcript_end = as.integer(transcript_end)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(start = min(transcript_start), end = min(transcript_end)) %>%
  dplyr::mutate(strand = 1) %>%
  dplyr::ungroup() %>%
  as.data.frame()
rownames(intron_metadata) = intron_metadata$transcript_id

#Calculate abundance ratios
gene_name_map = dplyr::select(intron_metadata, gene_id, transcript_id)
intron_ratios = calculateTranscriptRatios(intron_counts, gene_name_map)

#Construct a SummarizedExperiment object
se = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = intron_counts, tpm_ratios = intron_ratios), 
  colData = colData(ensembl_quants), 
  rowData = intron_metadata)
saveRDS(se, "results/acLDL/acLDL_leafcutter_counts.rds")
