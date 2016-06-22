library("DESeq2")
library("dplyr")
library("devtools")
library("gplots")
library("ggplot2")
library("tidyr")
load_all("../seqUtils/")
library("pheatmap")

#Load processed expression data from disk
expression_data = readRDS("results/acLDL/acLDL_combined_expression_data.rds")
weak_responders = c("xugn","giuo","oarz","fikt","kuxp","oefg","coio","nusw","hiaf","voas","cicb","auim","eiwy")
very_weak_responders = c("voas","coio","giuo", "oefg","oarz", "hiaf","kuxp","piun", "xugn","cicb","fikt", "nusw")

expression_data$sample_metadata = dplyr::mutate(expression_data$sample_metadata, response_group = ifelse(donor %in% weak_responders, "weak", "strong"))

#Import PEER residuals
peer_residuals = readRDS("results/acLDL/PEER/output/PEER_residuals.rds")

#Find expressed genes
mean_expression = seqUtils::calculateMean(expression_data$cqn, as.data.frame(expression_data$sample_metadata), "condition_name")
expressed_genes = names(which(apply(mean_expression, 1, max) > 0.5))
expressed_data = extractGenesFromExpressionList(expression_data, expressed_genes)
expressed_data_2 = extractConditionFromExpressionList(c("Ctrl","AcLDL"), expressed_data)

#Create purity df
purity = dplyr::select(expressed_data$sample_metadata, sample_id, mean_purity) %>% as.data.frame()
rownames(purity) = purity$sample_id
purity = dplyr::select(purity, -sample_id)

#Make heatmap of gene expression
cor_matrix = cor(expressed_data_2$cqn, method = "spearman")
pheatmap(cor_matrix, annotation_col = purity)
pheatmap(cor_matrix, annotation_col = purity, filename = "results/acLDL/DE/expressed_genes_heatmap.pdf", width = 16, height = 14,border_color = NA)

#Look at PEER residuals
cor_matrix = cor(peer_residuals, method = "pearson")
pheatmap(cor_matrix, filename ="results/acLDL/DE/expressed_genes.PEER_residuals.pdf", width = 16, height = 14,border_color = NA)


#Perform PCA
pca_list = performPCA(expressed_data$cqn[de_results_filtered$gene_id,], expressed_data$sample_metadata)
pca_plot = ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition, label = sample_id)) + 
  geom_point() +
  geom_text()
ggsave("results/acLDL/DE/expressed_genes_PCA.pdf", pca_plot, width = 12, height = 10)

#Perform PCA on AcLDL condition only
acLDL_samples = dplyr::filter(expressed_data$sample_metadata, condition_name == "AcLDL")
pca_list = performPCA(expressed_data$cqn[de_results_filtered$gene_id,acLDL_samples$sample_id], acLDL_samples)
pca_plot = ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition, label = sample_id)) + 
  geom_point() +
  geom_text()
#Extract gene metadata
gene_data = dplyr::select(expressed_data_2$gene_metadata, gene_id, gene_name, gene_biotype)

### AcLDL vs Ctrl DE ###
#Construct design matrix
#ctrl_vs_acldl_design = dplyr::filter(expressed_data_2$sample_metadata, condition_name %in% c("AcLDL","Ctrl"), !(donor %in% c("piun")))
#ctrl_vs_acldl_design = dplyr::filter(expressed_data_2$sample_metadata, condition_name %in% c("AcLDL","Ctrl"))
ctrl_vs_acldl_design = dplyr::filter(expressed_data_2$sample_metadata, condition_name %in% c("AcLDL","Ctrl")) %>%
  dplyr::filter(!(donor %in% very_weak_responders))
rownames(ctrl_vs_acldl_design) = ctrl_vs_acldl_design$sample_id
filtered_counts = expressed_data$counts[,ctrl_vs_acldl_design$sample_id]

#Peform differential expression analysis together with response groups
#dds = DESeqDataSetFromMatrix(filtered_counts, ctrl_vs_acldl_design, ~condition + response_group + condition*response_group)
#dds = DESeq(dds)
#res = results(dds) #Extract DE results

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
  arrange(-log2FoldChange)
write.table(de_results_filtered, "results/acLDL/DE/ctrl_vs_acldl_DE_genes_filtered.txt", sep= "\t", quote = FALSE, row.names = FALSE)

#Make heatmap with DE genes only (1.5 fold)
cor_matrix = cor(expressed_data$cqn[de_results_filtered$gene_id, ctrl_vs_acldl_design$sample_id], method = "pearson")
pheatmap(cor_matrix, annotation_col = purity, filename = "results/acLDL/DE/58_sample_1.5fold_heatmap.pdf", width = 15, height = 12, border_color = NA)

pca_list = performPCA(expressed_data$cqn[de_results_filtered$gene_id, ctrl_vs_acldl_design$sample_id], ctrl_vs_acldl_design)
ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition, label = sample_id)) + 
  geom_point() +
  geom_text()

pca_list = performPCA(expressed_data$cqn[de_results_filtered$gene_id, dplyr::filter(ctrl_vs_acldl_design, condition_name == "AcLDL")$sample_id], ctrl_vs_acldl_design)
ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition, label = sample_id)) + 
  geom_point() +
  geom_text()

cor_matrix = cor(expressed_data$cqn[de_results_filtered$gene_id,ctrl_vs_acldl_design_full$sample_id], method = "spearman")
pheatmap(cor_matrix, annotation_col = purity)

pheatmap(cor_matrix, annotation_col = purity, filename = "results/acLDL/DE/DE_genes_1,5fold_heatmap.pdf", width = 16, height = 14,border_color = NA)

ctrl_vs_acldl_design = dplyr::filter(expressed_data_2$sample_metadata, condition_name %in% c("AcLDL","Ctrl")) %>%
  dplyr::filter(!(donor %in% very_weak_responders))
ctrl_vs_acldl_design_full = dplyr::filter(expressed_data_2$sample_metadata, condition_name %in% c("AcLDL","Ctrl"))

#Look at PEER residuals
cor_matrix = cor(peer_residuals[intersect(rownames(peer_residuals),de_results_filtered$gene_id), ], method = "pearson")
pheatmap(cor_matrix, filename ="results/acLDL/DE/DE_genes_1,5fold_heatmap.PEER_residuals.pdf", width = 16, height = 14,border_color = NA)

#2-fold change in DE
de_results_filtered2 = dplyr::filter(de_results, padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(-log2FoldChange)
cor_matrix = cor(expressed_data$cqn[de_results_filtered2$gene_id,ctrl_vs_acldl_design$sample_id], method = "pearson")
pheatmap(cor_matrix, annotation_col = purity)
pheatmap(cor_matrix, annotation_col = purity, filename = "results/acLDL/DE/58_sample_2fold_heatmap.pdf", width = 15, height = 12, border_color = NA)

cor_matrix = cor(expressed_data$cqn[de_results_filtered2$gene_id,ctrl_vs_acldl_design_full$sample_id], method = "spearman")
pheatmap(cor_matrix, annotation_col = purity)


pheatmap(cor_matrix, annotation_col = purity, filename = "results/acLDL/DE/DE_genes_2fold_heatmap.pdf", width = 16, height = 14, border_color = NA)

#Look at PEER residuals
cor_matrix = cor(peer_residuals[intersect(rownames(peer_residuals),de_results_filtered2$gene_id), ], method = "pearson")
pheatmap(cor_matrix,filename = "results/acLDL/DE/DE_genes_2fold_heatmap.PEER_residuals.pdf",annotation_col = purity,width = 16, height = 14, border_color = NA)

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

#Export counts matrix to disk
write.table(expression_data$counts, "results/acLDL/DE/acLDL_counts_matrix.txt", sep = "\t", quote = FALSE)
write.table(expression_data$cqn, "results/acLDL/DE/acLDL_cqn_matrix.txt", sep = "\t", quote = FALSE)


#Explore the effect of weak responders
acldl_meta = dplyr::filter(expression_data$sample_metadata, condition_name == "AcLDL")
purity = ggplot(acldl_meta, aes(x = response_group, y = mean_purity_filtered)) + geom_boxplot() + geom_jitter(position = position_jitter(width = 0.1))
ggsave(filename = "results/acLDL/DE/response_group_purity.pdf", plot = purity, width = 6, height = 6)

stimulation_date = ggplot(acldl_meta, aes(y = response_group, x = factor(acLDL_date))) + geom_jitter(position=position_jitter(height=0.2, width=0))
ggsave(filename = "results/acLDL/DE/response_group_date.pdf", plot = stimulation_date, width = 9, height = 6)

rna_conc = ggplot(acldl_meta, aes(x = response_group, y = ng_ul_mean)) + geom_boxplot() + geom_jitter(position = position_jitter(width = 0.1))
ggsave(filename = "results/acLDL/DE/response_group_rna_conc.pdf", plot = rna_conc, width = 6, height = 6)


