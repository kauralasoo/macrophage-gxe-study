library("DESeq2")
library("dplyr")
library("devtools")
library("gplots")
library("ggplot2")
library("tidyr")
load_all("../seqUtils/")

#Load processed expression data from disk
expression_data = readRDS("results/acLDL/acLDL_combined_expression_data.rds")

#Find expressed genes
mean_expression = seqUtils::calculateMean(expression_data$cqn, as.data.frame(expression_data$sample_metadata), "condition_name")
expressed_genes = names(which(apply(mean_expression, 1, max) > 0.5))
expressed_data = extractGenesFromExpressionList(expression_data, expressed_genes)

#Make heatmap of gene expression
cor_matrix = cor(expressed_data$cqn, method = "spearman")
heatmap.2(cor_matrix, margins = c(12,12))

#Perform PCA
pca_list = performPCA(expressed_data$cqn, expressed_data$sample_metadata)
ggplot(pca_list$pca_matrix, aes(x = PC4, y = PC5, color = condition, label = sample_id)) + 
  geom_point() +
  geom_text()

#Peform differential expression analysis
dds = DESeqDataSetFromMatrix(filtered_counts, filtered_design, ~condition)
dds = DESeq(dds)
res = results(dds) #Extract DE results
de_results = res %>%
  as.data.frame() %>% #Convert to data.frame
  dplyr::mutate(gene_id = rownames(res)) %>% #Add gene id as a column
  dplyr::left_join(gene_data, by = "gene_id") %>% #AD
  tbl_df %>% #Convert to a convenient dplyr class
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 1) %>% 
  dplyr::arrange(desc(log2FoldChange))

#Don't apply any fold-change filtering
de_results_unfiltered = res %>%
  as.data.frame() %>% #Convert to data.frame
  dplyr::mutate(gene_id = rownames(res)) %>% #Add gene id as a column
  dplyr::left_join(gene_data, by = "gene_id") %>% #AD
  tbl_df()

#Save results to disk
write.table(de_results, "results/acLDL/DE/DE_genes.txt", sep= "\t", quote = FALSE, row.names = FALSE)
write.table(de_results_unfiltered, "results/acLDL/DE/DE_genes_unfiltered.txt", sep= "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_exprs_cqn, "results/acLDL/expression_cqn_normalised.txt", sep = "\t", quote = FALSE)
