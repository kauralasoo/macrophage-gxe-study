library("readr")
library("devtools")
library("dplyr")
library("tidyr")
library("SummarizedExperiment")
load_all("../seqUtils/")

#Import read counts
combined_expression_data = readRDS("results/SL1344/combined_expression_data.rds")
sample_ids = c("chrom", combined_expression_data$sample_metadata$sample_id)

#Construct sample metadata df
sample_metadata = combined_expression_data$sample_metadata %>% as.data.frame()
rownames(sample_metadata) = sample_metadata$sample_id

#Import intron and cluster counts from leafcutter
intron_counts = readr::read_delim("results/SL1344/leafcutter/leafcutter.intron_counts.txt.gz", 
                                  delim = " ", col_names = TRUE)[,sample_ids]
cluster_counts = readr::read_delim("results/SL1344/leafcutter/leafcutter.cluster_counts.txt.gz", 
                                   delim = " ", col_names = TRUE)[,sample_ids]

#Convert to matrices
intron_matrix = dplyr::select(intron_counts, -chrom) %>% as.matrix()
rownames(intron_matrix) = intron_counts$chrom

cluster_matrix = dplyr::select(cluster_counts, -chrom) %>% as.matrix()
rownames(cluster_matrix) = cluster_counts$chrom

#Calculate ratios
intron_ratios = intron_matrix/cluster_matrix

#Extract intron coords from counts
intron_coords = dplyr::select(intron_counts, chrom) %>%
  tidyr::separate(chrom, c("chr","intron_start","intron_end","cluster_id"), sep = ":", remove = FALSE) %>%
  dplyr::mutate(intron_start = as.integer(intron_start), intron_end = as.integer(intron_end)) %>%
  dplyr::rename(transcript_id = chrom) %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::mutate(cluster_start = min(intron_start), cluster_end = max(intron_end), cluster_size = length(cluster_id)) %>%
  dplyr::ungroup() %>%
  as.data.frame()
rownames(intron_coords) = intron_coords$transcript_id

se = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = intron_matrix, count_ratios = intron_ratios), 
  colData = sample_metadata, 
  rowData = intron_coords)

saveRDS(se,"results/SL1344/leafCutter_summarized_experiment.rds")



