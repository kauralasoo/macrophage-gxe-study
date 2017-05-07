library("dplyr")
library("devtools")
library("GenomicRanges")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../seqUtils/")

#Import Ensembl quant results for reference
ensembl_quants = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
sample_names = colnames(ensembl_quants)

#Import Ensembl transcript database
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names = TRUE)

#Import transcript metadata
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")
transcript_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.transcript_data.rds")) %>%
  dplyr::filter(gene_biotype %in% valid_gene_biotypes, transcript_gencode_basic == "GENCODE basic", 
                chromosome_name %in% valid_chromosomes)
filtered_exons = exons[transcript_data$ensembl_transcript_id]

#Import intron and cluster counts from leafcutter
intron_counts = read.table("processed/acLDL/leafcutter/leafcutter_perind_numers.counts.gz")[,sample_names]

#Extract intron coords from counts
intron_metadata = leafcutterConstructMeta(rownames(intron_counts))

#Convert leafcutter clusters into a Granges list for visualisation
granges = leafcutterMetaToGrangesList(intron_metadata)
leafcutter_grl = GRangesList(granges)
saveRDS(leafcutter_grl, "results/acLDL/acLDL_leafcutter.GRangesList.rds")

#Find overlas between two sets of exons
overlaps = findOverlaps(leafcutter_grl, filtered_exons, ignore.strand = TRUE)
matches = data_frame(transcript_id = names(leafcutter_grl[queryHits(overlaps)]),
                     ensembl_transcript_id = names(filtered_exons[subjectHits(overlaps)]))

#Add gene names to overlaps
filtered_matches = dplyr::left_join(matches, dplyr::select(intron_metadata, transcript_id, gene_id), by = "transcript_id") %>%
  dplyr::left_join(dplyr::select(transcript_data, ensembl_transcript_id, ensembl_gene_id, external_gene_name), by = "ensembl_transcript_id")

#Collapse to cluster level
gene_matches = dplyr::transmute(filtered_matches, gene_id, ensembl_gene_id, gene_name = external_gene_name) %>%
  unique()
grouped_matches = dplyr::group_by(gene_matches, gene_id) %>%
  dplyr::summarise(ensembl_gene_id = paste(ensembl_gene_id, collapse =";"),
                   gene_name = paste(gene_name, collapse = ";"), 
                   overlap_count = length(gene_id)) %>%
  dplyr::ungroup()

#Add gene strand to grouped matches
strand_info = dplyr::transmute(transcript_data, ensembl_gene_id, ensembl_gene_strand = ifelse(strand == 1, "+", "-")) %>% unique()
final_matches = dplyr::left_join(grouped_matches, strand_info, by = "ensembl_gene_id")

#Compile final intron metadata
intron_meta_final = dplyr::left_join(intron_metadata, final_matches, by = "gene_id") %>%
  as.data.frame()
rownames(intron_meta_final) = intron_metadata$transcript_id

#Calculate abundance ratios
gene_name_map = dplyr::select(intron_metadata, gene_id, transcript_id)
intron_ratios = calculateTranscriptRatios(intron_counts, gene_name_map)

#Construct a SummarizedExperiment object
se = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = intron_counts, tpm_ratios = intron_ratios), 
  colData = colData(ensembl_quants), 
  rowData = intron_meta_final)
saveRDS(se, "results/acLDL/acLDL_leafcutter_counts.rds")

#Make a new Granges list based on overlapping gene strand

#Make Granges list with the correct strand
leafcutter_ranges = dplyr::transmute(tbl_df(intron_meta_final), transcript_id, chr, transcript_start, 
                                     transcript_end, strand = ifelse(!is.na(ensembl_gene_strand), ensembl_gene_strand, strand)) %>% 
  leafcutterMetaToGrangesList() %>%
  GRangesList()
saveRDS(leafcutter_ranges, "results/acLDL/acLDL_leafcutter.GRangesList.rds")
