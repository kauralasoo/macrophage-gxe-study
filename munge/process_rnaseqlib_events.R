library("dplyr")
library("tidyr")
library("devtools")
library("rtracklayer")
library("wiggleplotr")
library("rtracklayer")
load_all("../reviseAnnotations/")

#Import events
RI_events = import.gff3("../../annotations/rnaseqlib/IRF5_gff/commonshortest/SE.hg38.gff3")
RI_list = revisedGffToGrangesList(RI_events)
plotTranscripts(RI_list, rescale_introns = FALSE)

#Construct event metadata
transcript_data = data_frame(transcript_id = names(RI_list)) %>%
  dplyr::mutate(split_id = strsplit(transcript_id, "\\:\\+")) %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::mutate(split_id_length = length(unlist(split_id))) %>%
  dplyr::mutate(first_coord = unlist(split_id)[1]) %>%
  dplyr::mutate(second_coord = unlist(split_id)[split_id_length-1]) %>%
  dplyr::ungroup() %>%
  tidyr::separate(first_coord, c("chr_name", "transcript_start", "f_end"), sep = ":") %>%
  tidyr::separate(chr_name, c("none", "chr"), sep = "chr") %>%
  tidyr::separate(second_coord, c("s_chr", "s_start", "transcript_end"), sep = ":") %>%
  dplyr::select(transcript_id, chr, transcript_start, transcript_end) %>%
  dplyr::mutate(event_id = paste("evt", chr, transcript_start, transcript_end, sep = "_")) %>%
  dplyr::mutate(transcript_start = as.integer(transcript_start), transcript_end = as.integer(transcript_end))
  
#Extract event data from transcript data
event_data = dplyr::select(transcript_data, event_id, chr, transcript_start, transcript_end) %>%
  unique()

event_ranges = dplyr::transmute(event_data, seqnames = chr, start = transcript_start, end = transcript_end, strand = "+") %>% 
  dataFrameToGRanges()
event_olaps = findOverlaps(event_ranges)

#Construct a graph
edges = event_olaps %>% as.data.frame() %>% 
  dplyr::filter(queryHits != subjectHits) %>% as.matrix() %>% t() %>% as.vector()
edges = c(edges, 2,3)
event_overlap_graph = igraph::graph(edges, directed = FALSE)

cl = igraph::clusters(event_overlap_graph)
grpahs = igraph::decompose.graph(event_overlap_graph)

igraph::max_cliques(event_overlap_graph)

dist_matrix = igraph::distances(event_overlap_graph)
rownames(dist_matrix) = event_data$event_id
colnames(dist_matrix) = event_data$event_id

