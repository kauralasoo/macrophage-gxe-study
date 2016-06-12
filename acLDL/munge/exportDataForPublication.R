library("dplyr")

acldl_data = readRDS("results/acLDL/acLDL_combined_expression_data.rds")

#Export sample metadata
write.table(acldl_data$sample_metadata, "results/acLDL/datasets/acldl_sample_metadata.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#Export open access samples
oa_samples = dplyr::filter(acldl_data$sample_metadata, open_access == 1) %>% 
  dplyr::select(sample_id, genotype_id)
write.table(oa_samples, "results/acLDL/datasets/acldl_open_access_samples.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#Extract OA column names from the ASE data
columns = read.table("results/acLDL/datasets/acldl_column_names.txt") %>% t() %>% as.vector()
column_df = data_frame(column_name = columns)
column_df$column_number = 1:nrow(column_df)

#extract columns
selected_columns = c(column_df$column_name[1:5], oa_samples$sample_id)
filtered_columns = dplyr::filter(column_df, column_name %in% selected_columns)

