library("dplyr")
library("SummarizedExperiment")
library("devtools")
load_all("../seqUtils/")

#Import Ensembl quant results for reference
ensembl_quants = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
sample_names = colnames(ensembl_quants)

#Import Ensembl transcript database
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
transcript_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.transcript_data.rds")) %>%
  dplyr::select(ensembl_gene_id, ensembl_transcript_id)

#Import intron and cluster counts from leafcutter
intron_counts = read.table("processed/acLDL/leafcutter/leafcutter_perind_numers.counts.gz")[,sample_names]

#Extract intron coords from counts
intron_metadata = dplyr::data_frame(transcript_id = rownames(intron_counts)) %>%
  tidyr::separate(transcript_id, c("chr","transcript_start","transcript_end","gene_id"), sep = ":", remove = FALSE) %>%
  dplyr::mutate(transcript_start = as.integer(transcript_start), transcript_end = as.integer(transcript_end)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(start = min(transcript_start), end = min(transcript_end)) %>%
  dplyr::mutate(strand = 1) %>%
  dplyr::ungroup()

#Find potential overlaps with Ensembl genes
leafcutter_clusters = 
  dplyr::transmute(intron_metadata, transcript_id, gene_id, 
                   seqnames = chr, start = transcript_start, 
                   end = transcript_end, strand = ifelse(strand == 1, "+", "-")) %>% 
  dataFrameToGRanges()
potential_overlaps = findOverlaps(leafcutter_clusters, exons, ignore.strand=TRUE)
matches = data_frame(ensembl_transcript_id = names(exons[subjectHits(potential_overlaps)]), 
                     gene_id = elementMetadata(leafcutter_clusters[queryHits(potential_overlaps)])$gene_id) %>%
  dplyr::left_join(transcript_data, by = "ensembl_transcript_id") %>%
  dplyr::select(gene_id, ensembl_gene_id) %>%
  unique()
grouped_matches = dplyr::group_by(matches, gene_id) %>%
  dplyr::summarise(ensembl_gene_id = paste(ensembl_gene_id, collapse =";")) %>%
  dplyr::ungroup()

#Compile final intron metadata
intron_metadata = dplyr::left_join(intron_metadata, grouped_matches, by = "gene_id") %>%
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
