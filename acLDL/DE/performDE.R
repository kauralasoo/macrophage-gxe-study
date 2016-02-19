library("DESeq2")
library("dplyr")
library("devtools")
library("gplots")
library("ggplot2")
library("tidyr")
load_all("../seqUtils/")

#Load processed expression data from disk
expression_data = readRDS("results/acLDL/acLDL_combined_expression_data.rds")
strong_responder_donors = scan("macrophage-gxe-study/data/sample_lists/acLDL/acLDL_strong_responders.txt", what = "character")
expression_data$sample_metadata =  dplyr::mutate(expression_data$sample_metadata, response_group = ifelse(donor %in% strong_responder_donors, "group_1", "group_2"))

#Find expressed genes
mean_expression = seqUtils::calculateMean(expression_data$cqn, as.data.frame(expression_data$sample_metadata), "condition_name")
expressed_genes = names(which(apply(mean_expression, 1, max) > 0.5))
expressed_data = extractGenesFromExpressionList(expression_data, expressed_genes)
expressed_data_2 = extractConditionFromExpressionList(c("Ctrl","AcLDL"), expressed_data)

#Import list of strong responders
strong_responders_meta = dplyr::filter(expressed_data$sample_metadata, donor %in% strong_responder_donors)
weak_responders_meta = dplyr::filter(expressed_data$sample_metadata, !(donor %in% strong_responder_donors))

#Make heatmap of gene expression
cor_matrix = cor(expressed_data$cqn[,strong_responders_meta$sample_id], method = "spearman")
heatmap.2(cor_matrix, margins = c(12,12))

#Perform PCA
pca_list = performPCA(expressed_data$cqn, expressed_data$sample_metadata)
ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition, label = sample_id)) + 
  geom_point() +
  geom_text()

#Extract gene metadata
gene_data = dplyr::select(expressed_data$gene_metadata, gene_id, gene_name, gene_biotype)

### AcLDL vs Ctrl DE ###
#Construct design matrix
ctrl_vs_acldl_design = dplyr::filter(expressed_data$sample_metadata, condition_name %in% c("AcLDL","Ctrl"))
rownames(ctrl_vs_acldl_design) = ctrl_vs_acldl_design$sample_id
filtered_counts = expressed_data$counts[,ctrl_vs_acldl_design$sample_id]

#Peform differential expression analysis together with response groups
dds = DESeqDataSetFromMatrix(filtered_counts, ctrl_vs_acldl_design, ~condition + response_group + condition*response_group)
dds = DESeq(dds)
res = results(dds) #Extract DE results

#Peform differential expression analysis
dds = DESeqDataSetFromMatrix(filtered_counts, ctrl_vs_acldl_design, ~condition)
dds = DESeq(dds)
res = results(dds) #Extract DE results
de_results = res %>%
  as.data.frame() %>% #Convert to data.frame
  dplyr::mutate(gene_id = rownames(res)) %>% #Add gene id as a column
  dplyr::left_join(gene_data, by = "gene_id") %>% #AD
  tbl_df() %>% #Convert to a convenient dplyr class
  dplyr::arrange(padj) %>%
  dplyr::select(gene_id, gene_name, everything())
write.table(de_results, "results/acLDL/DE/ctrl_vs_acldl_DE_genes.txt", sep= "\t", quote = FALSE, row.names = FALSE)

#Filter by p-value and fold-change
de_results_filtered = dplyr::filter(de_results, padj < 0.05, abs(log2FoldChange) > 0.58) %>%
  arrange(desc(log2FoldChange))
write.table(de_results_filtered, "results/acLDL/DE/ctrl_vs_acldl_DE_genes_filtered.txt", sep= "\t", quote = FALSE, row.names = FALSE)

#Make heatmap with DE genes only
de_results_filtered$gene_id

cor_matrix = cor(expressed_data$cqn[de_results_filtered$gene_id,weak_responders_meta$sample_id ], method = "pearson")
heatmap.2(cor_matrix, margins = c(12,12))

#Test for interaction between LAL and AcLDL
lal_design = dplyr::semi_join(expressed_data$sample_metadata, 
                 dplyr::filter(expressed_data$sample_metadata, condition_name %in% c("LAL")), by = "donor") %>%
  dplyr::select(sample_id, donor, condition_name) %>%
  dplyr::mutate(AcLDL = ifelse(condition_name %in% c("AcLDL","LAL_AcLDL"), "yes", "no")) %>%
  dplyr::mutate(LAL = ifelse(condition_name %in% c("LAL","LAL_AcLDL"), "yes","no"))
rownames(lal_design) = lal_design$sample_id

lal_counts = expressed_data$counts[,lal_design$sample_id]

#Test DE
dds = DESeqDataSetFromMatrix(lal_counts, lal_design, ~AcLDL + LAL + AcLDL*LAL)
dds = DESeq(dds)
results(dds, name = "AcLDLyes.LALyes") %>% as.data.frame() %>% dplyr::filter(padj < 0.1)
results(dds, name = "AcLDL_yes_vs_no") %>% as.data.frame() %>% dplyr::filter(padj < 0.05)
results(dds, name = "LAL_yes_vs_no") %>% as.data.frame() %>% dplyr::filter(padj < 0.05)

