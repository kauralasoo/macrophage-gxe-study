library("devtools")
library("cqn")
load_all("macrophage-gxe-study/seqUtils/")

#Import raw read counts 
data = read.table("results/SL1344/combined_counts.txt", stringsAsFactors = FALSE)
#There seems to be a sample swap between coio_C and coio_D, fix that
indexes = c(which(colnames(data) == "coio_C"), which(colnames(data) == "coio_D"))
colnames(data)[indexes] = c("coio_D", "coio_C")

#Load transcript metadata
filtered_metadata = readRDS("annotations/biomart_transcripts.filtered.rds")
gene_data = dplyr::select(filtered_metadata, gene_id, gene_name, gene_biotype, percentage_gc_content) %>% unique()
length_df = dplyr::select(data, gene_id, length)
gene_metadata = dplyr::left_join(gene_data, length_df, by = "gene_id") %>% tbl_df()

#Process counts
filtered_data = dplyr::filter(data, gene_id %in% gene_data$gene_id)
counts = dplyr::select(filtered_data, -gene_id, -length)
rownames(counts) = filtered_data$gene_id
counts = counts[gene_metadata$gene_id,] #Reoder counts

#Construct a design matrix from the sample names
design_matrix = constructDesignMatrix_SL1344(sample_ids = colnames(counts))

#cqn_normalize_data
exprs_cqn = calculateCQN(counts, gene_metadata)

results_list = list(
  exprs_counts = counts,
  exprs_cqn = exprs_cqn,
  design = design_matrix,
  gene_metadata = gene_metadata)

saveRDS(results_list, "results/SL1344/combined_expression_data.rds")
