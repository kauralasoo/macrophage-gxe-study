#Export cleaned up and processed data for sharing

#Import ATAC data
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
atac_data_no_covariates = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Export peak accessibility data
#Counts
peak_access_counts <- gzfile("results/ATAC/clean_data/peak_accessibility_counts.txt.gz", "w")
write.table(atac_data$counts, peak_access_counts, sep = "\t", quote = FALSE)
close(peak_access_counts)

#CQN
peak_access_cqn <- gzfile("results/ATAC/clean_data/peak_accessibility_cqn.txt.gz", "w")
write.table(atac_data$cqn, peak_access_cqn, sep = "\t", quote = FALSE)
close(peak_access_cqn)

#Export peak metadata
write.table(atac_data$gene_metadata, "results/ATAC/clean_data/peak_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Export sample metadata
write.table(atac_data_no_covariates$sample_metadata, "results/ATAC/clean_data/sample_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Export sample metadata before filtering (to highlight QC criteria)
atac_qc_metadata = readRDS("macrophage-gxe-study/data/chromatin/ATAC/compiled_atac_metadata.rds") %>%
  dplyr::transmute(sample_id, donor, condition = condition_name, assigned, assigned_frac = round(assigned_frac,2), 
                   mt_frac = round(mt_frac,2), duplication = round(percent_duplication,2), peak_count, length_ratio = round(short_long_ratio,2)) %>%
  dplyr::mutate(qtl_mapping = ifelse(sample_id %in% c("fikt_B_ATAC", "qaqx_A_ATAC", "qaqx_B_ATAC", "uaqe_A_ATAC", "uaqe_B_ATAC"), "no","yes")) %>%
  dplyr::mutate(diff_access = ifelse(donor %in% 
    c("qolg","vass","kuxp","cicb", "febc","eiwy","oapg","nukw","hayt","bima","pamv","guss","eipl","iill","podx","pelm"),"yes","no")) %>%
  dplyr::select(-donor) %>%
  dplyr::arrange(sample_id)
write.table(atac_qc_metadata, "results/ATAC/clean_data/atac_qc_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)
