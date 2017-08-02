library("rtracklayer")

#Define helper functions
findValidPairs <- function(master_id, meta_df, bin_df){
  
  #Shuffle metadata
  meta_df = meta_df[sample(nrow(meta_df),nrow(meta_df)),]
  
  #Extract master record
  master_record = dplyr::filter(meta_df, gene_id == master_id)
  
  #Find potential dependent peaks
  potential_dependents = dplyr::semi_join(meta_df, master_record, by = "chr") %>%
    dplyr::transmute(master_id = as.character(master_record$gene_id), dependent_id = gene_id)
  
  #Calculate distances between master and dependent peaks
  dist_table = calculatePeakDistance(potential_dependents, meta_df) %>%
    dplyr::mutate(distance = abs(distance)) %>%
    dplyr::filter(distance > 0) %>%
    dplyr::filter(distance < 55000)
  
  valid_pairs = dplyr::mutate(dist_table, bin_id = ((distance %/% 2500)+1) * 2500) %>% 
    dplyr::left_join(bin_df, by = "bin_id") %>% 
    dplyr::mutate(prob = runif(length(fraction))) %>% 
    dplyr::filter(prob < fraction)
  
  return(valid_pairs)
}

#Import atac data
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Import Master-dependent pairs
result_list = readRDS("results/ATAC/QTLs/qtl_peak_type_assignment.rds")
master_dependent_pairs = result_list$dependents$unique_masters

#Calculate real distances
dependent_distances = calculatePeakDistance(dplyr::select(master_dependent_pairs, master_id, dependent_id), atac_list$gene_metadata)
histogram = hist(abs(dependent_distances$distance), breaks = seq(0,55000,2500))
bins = histogram$breaks[-1]
bin_frequencies = histogram$counts/sum(histogram$counts)
bin_df = data_frame(bin_id = bins, fraction = bin_frequencies)

#Construct random pairs
meta_df = tbl_df(atac_list$gene_metadata)

#Construct random pairs
sample1 = meta_df[sample(nrow(meta_df),nrow(meta_df)),]
master_ids = sample1$gene_id[1:5000] %>% idVectorToList()
valid_pairs = purrr::map_df(master_ids, ~findValidPairs(.,meta_df, bin_df))
saveRDS(valid_pairs, "results/ATAC/ATAC_random_pairs.rds")

#Construct pairs
peak_coords_gr = dplyr::bind_rows(dplyr::transmute(valid_pairs, peak_id = master_id, seqnames = master_chr, start = master_centre), 
                 dplyr::transmute(valid_pairs, peak_id = dependent_id, seqnames = dependent_chr, start = dependent_centre)) %>%
  unique() %>%
  dplyr::mutate(end = start, strand = "*") %>%
  dplyr::mutate(seqnames = paste0("chr",seqnames)) %>%
  dataFrameToGRanges()

#Lift over coordinates using rtracklayer::liftOver
chain = rtracklayer::import.chain("macrophage-gxe-study/data/liftOver_genotypes/hg38ToHg19.over.chain")
new_coords = rtracklayer::liftOver(peak_coords_gr, chain) %>% as.list()
new_coords_df = purrr::map_df(new_coords, ~as.data.frame(.)) %>% tbl_df()
new_coords_dff =tidyr::separate(new_coords_df, seqnames, c("none", "chr"), sep = "chr") %>%
  dplyr::select(peak_id, chr, start)

#Add coords back
random_pairs = dplyr::select(valid_pairs, master_id, dependent_id) %>% 
  dplyr::left_join(new_coords_dff, by = c("master_id" = "peak_id")) %>% 
  dplyr::rename(master_chr = chr, master_midpoint = start) %>%
  dplyr::left_join(new_coords_dff, by = c("dependent_id" = "peak_id")) %>% 
  dplyr::rename(master_chr = chr, dependent_midpoint = start)
write.table(random_pairs, "results/ATAC/TAD_enrichment/random_multipeak_pairs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#Convert true pairs to GRCh37 coordinates
peak_coords_gr = dplyr::bind_rows(dplyr::transmute(dependent_distances, peak_id = master_id, seqnames = master_chr, start = master_centre), 
                                  dplyr::transmute(dependent_distances, peak_id = dependent_id, seqnames = dependent_chr, start = dependent_centre)) %>%
  unique() %>%
  dplyr::mutate(end = start, strand = "*") %>%
  dplyr::mutate(seqnames = paste0("chr",seqnames)) %>%
  dataFrameToGRanges()

#Lift over coordinates using rtracklayer::liftOver
chain = rtracklayer::import.chain("macrophage-gxe-study/data/liftOver_genotypes/hg38ToHg19.over.chain")
new_coords = rtracklayer::liftOver(peak_coords_gr, chain) %>% as.list()
new_coords_df = purrr::map_df(new_coords, ~as.data.frame(.)) %>% tbl_df()
new_coords_dff =tidyr::separate(new_coords_df, seqnames, c("none", "chr"), sep = "chr") %>%
  dplyr::select(peak_id, chr, start)

#Add coords back
true_pairs = dplyr::select(dependent_distances, master_id, dependent_id) %>% 
  dplyr::left_join(new_coords_dff, by = c("master_id" = "peak_id")) %>% 
  dplyr::rename(master_chr = chr, master_midpoint = start) %>%
  dplyr::left_join(new_coords_dff, by = c("dependent_id" = "peak_id")) %>% 
  dplyr::rename(master_chr = chr, dependent_midpoint = start)
write.table(true_pairs, "results/ATAC/TAD_enrichment/true_multipeak_pairs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#Make a plot of true distances
dependent_distance_plot = ggplot(dependent_distances, aes(x = abs(distance)/1000)) + geom_histogram(binwidth = 2) + theme_light() +
  xlab("Distance from master region (kb)") +
  ylab("Dependent region count")
ggsave("figures/main_figures/caQTL_master_dependent_peak_distance.pdf", plot = dependent_distance_plot, width = 3, height = 3)


#Import TAD overlap results from Natsuhiko
true_tad = read.table("results/ATAC/TAD_enrichment/true_tad.txt", stringsAsFactors = FALSE) %>% tbl_df()
colnames(true_tad) = c("master_id","dependent_id","master_chr","master_midpoint","dependent_chr","dependent_midpoint", "tad_id")
true_tad = dplyr::mutate(true_tad, is_in_tad = ifelse(is.na(tad_id), FALSE, TRUE)) %>%
  dplyr::mutate(distance = abs(master_midpoint-dependent_midpoint))

rand_tad = read.table("results/ATAC/TAD_enrichment/rand_tad.txt", stringsAsFactors = FALSE) %>% tbl_df()
colnames(rand_tad) = c("master_id","dependent_id","master_chr","master_midpoint","dependent_chr","dependent_midpoint", "tad_id")
rand_tad = dplyr::mutate(rand_tad, is_in_tad = ifelse(is.na(tad_id), FALSE, TRUE)) %>%
  dplyr::mutate(distance = abs(master_midpoint-dependent_midpoint))

table(true_tad$is_in_tad)[2]/sum(table(true_tad$is_in_tad)) #75%
table(rand_tad$is_in_tad)[2]/sum(table(rand_tad$is_in_tad)) #74%

