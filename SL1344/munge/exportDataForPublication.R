library("dplyr")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data.rds")
combined_expression_data_cov = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Export sample metadata
write.table(combined_expression_data$sample_metadata, "results/SL1344/datasets/RNA/rna_sample_metadata.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#Export open access samples
oa_samples = dplyr::filter(combined_expression_data$sample_metadata, open_access == 1) %>% 
  dplyr::select(sample_id, genotype_id)
write.table(oa_samples, "results/SL1344/datasets/RNA/rna_open_access_samples.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#Extract OA column names from the ASE data
columns = read.table("results/SL1344/datasets/RNA/count_column_names.txt") %>% t() %>% as.vector()
column_df = data_frame(column_name = columns)
column_df$column_number = 1:nrow(column_df)

#extract columns
selected_columns = c(column_df$column_name[1:5], oa_samples$sample_id)
filtered_columns = dplyr::filter(column_df, column_name %in% selected_columns)

#Export read counts
write.table(combined_expression_data$counts, "results/SL1344/datasets/RNA/rna_gene_read_counts.txt", 
            sep = "\t", quote = FALSE)

#Export gene metadata
write.table(combined_expression_data$gene_metadata, "results/SL1344/datasets/RNA/rna_gene_metadata.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)


#Export lead QTL variants
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds") %>%
  purrr::map(~dplyr::select(.,-p_fdr))
write.table(rasqual_min_pvalues$naive, "results/SL1344/datasets/eQTLs/eQTL_naive_RASQUAL_lead_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rasqual_min_pvalues$IFNg, "results/SL1344/datasets/eQTLs/eQTL_IFNg_RASQUAL_lead_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rasqual_min_pvalues$SL1344, "results/SL1344/datasets/eQTLs/eQTL_SL1344_RASQUAL_lead_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rasqual_min_pvalues$IFNg_SL1344, "results/SL1344/datasets/eQTLs/eQTL_IFNg_SL1344_RASQUAL_lead_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds") %>%
  purrr::map(~dplyr::select(.,-p_fdr))


fastqtl_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")
write.table(fastqtl_min_pvalues$naive, "results/SL1344/datasets/eQTLs/eQTL_naive_FastQTL_lead_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fastqtl_min_pvalues$IFNg, "results/SL1344/datasets/eQTLs/eQTL_IFNg_FastQTL_lead_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fastqtl_min_pvalues$SL1344, "results/SL1344/datasets/eQTLs/eQTL_SL1344_FastQTL_lead_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fastqtl_min_pvalues$IFNg_SL1344, "results/SL1344/datasets/eQTLs/eQTL_IFNg_SL1344_FastQTL_lead_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

