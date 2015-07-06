library("DESeq2")
library("dplyr")
library("devtools")
library("gplots")
library("ggplot2")
library("tidyr")
load_all("macrophage-gxe-study/seqUtils/")

#Load processed expression data from disk
dataset = readRDS("results/acLDL/acLDL_combined_expression_data.rds")

#Load gene metadata
metadata = readRDS("annotations/Homo_sapiens.GRCh38.79.transcript_data.rds")
gene_data = dplyr::select(metadata, ensembl_gene_id, external_gene_name, gene_biotype) %>% unique() %>%
  dplyr::rename(gene_id = ensembl_gene_id, gene_name = external_gene_name)

#Remove outlier samples xegx and nusw
filtered_design = dplyr::filter(dataset$design, !(donor %in% c("nusw","xegx")))
rownames(filtered_design) = filtered_design$sample_id
filtered_exprs_cqn = dataset$exprs_cqn[,filtered_design$sample_id]
filtered_counts = dataset$exprs_counts[,filtered_design$sample_id]

#Find expressed genes
cqn_expressed_genes = names(which(rowMeans(filtered_exprs_cqn) > 0))
cqn_expressed = filtered_exprs_cqn[cqn_expressed_genes,]

#Make heatmap of gene expression
cor_matrix = cor(cqn_expressed, method = "spearman")
heatmap.2(cor_matrix, margins = c(12,12))


#Perform PCA
tpm_z = zScoreNormalize(tpm_expressed)
pca_list = performPCA(cqn_expressed, dataset$design)
ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition, label = sample_id)) + 
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

#Save results to disk
write.table(de_results, "results/acLDL/DE_genes.txt", sep= "\t", quote = FALSE, row.names = FALSE)
