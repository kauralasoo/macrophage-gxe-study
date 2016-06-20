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

