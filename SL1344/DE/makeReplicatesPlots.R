library("devtools")
library("cqn")
library("dplyr")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

replicates_data = readRDS("results/SL1344/replicates_combined_expression_data.rds")
other_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#Keep expressed genes
replicates_data = extractGenesFromExpressionList(replicates_data, rownames(other_data$cqn))

#Make some heatmaps
meta = as.data.frame(dplyr::select(replicates_data$sample_metadata, donor)) %>%
  dplyr::mutate(donor = ifelse(donor %in% c("nibo","fpdj"), "ffdj",donor))
rownames(meta) = replicates_data$sample_metadata$sample_id

pheatmap(cor(replicates_data$cqn, method = "spearman"), annotation_row = meta, 
         show_rownames = FALSE, show_colnames = TRUE)

#Filter by line variability
naive_meta = meta[meta$condition_name == "naive",]
pheatmap(cor(replicates_data$cqn[,rownames(naive_meta)], method = "spearman"), annotation_row = meta, 
         show_rownames = FALSE, show_colnames = TRUE)
pheatmap(cor(replicates_data$cqn[,rownames(naive_meta)], method = "spearman"), annotation_col = meta, 
         show_rownames = TRUE,show_colnames = FALSE, filename = "figures/supplementary/replicability_naive_heatmap.pdf", width = 5.5, height = 4)
