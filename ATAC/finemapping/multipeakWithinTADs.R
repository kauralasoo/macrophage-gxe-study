
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


dependent_distance_plot = ggplot(dependent_distances, aes(x = abs(distance)/1000)) + geom_histogram(binwidth = 2) + theme_light() +
  xlab("Distance from master region (kb)") +
  ylab("Dependent region count")
ggsave("figures/main_figures/caQTL_master_dependent_peak_distance.pdf", plot = dependent_distance_plot, width = 3, height = 3)

