library("DESeq2")
library("dplyr")
library("devtools")
library("tidyr")
library("ggplot2")
library("gplots")
library("cqn")
load_all("macrophage-gxe-study/seqUtils/")

#Load data from disk
data = read.table("results//SL1344//combined_counts.txt", stringsAsFactors = FALSE)
#There seems to be a sample swap between coio_C and coio_D, fix that
indexes = c(which(colnames(data) == "coio_C"), which(colnames(data) == "coio_D"))
colnames(data)[indexes] = c("coio_D", "coio_C")

#Load metadata
filtered_metadata = readRDS("annotations/biomart_transcripts.filtered.rds")
gene_data = dplyr::select(filtered_metadata, gene_id, gene_name, gene_biotype, percentage_gc_content) %>% unique()
length_df = dplyr::select(data, gene_id, length)
gene_metadata = dplyr::left_join(gene_data, length_df, by = "gene_id")

#Filter genes
filtered_data = dplyr::filter(data, gene_id %in% gene_data$gene_id)
counts = dplyr::select(filtered_data, -gene_id, -length)
rownames(counts) = filtered_data$gene_id
counts = counts[gene_metadata$gene_id,] #Reoder counts

#Construct a design matrix from the sample names
design_matrix = constructDesignMatrix_SL1344(sample_ids = colnames(counts))
design_filtered = dplyr::filter(design_matrix, donor %in% c("fpdj","vorx","zuta"))
rownames(design_filtered) = design_filtered$sample_id
counts = counts[,design_filtered$sample_id]

#CQN normalize the data
exprs_cqn = calculateCQN(counts, gene_metadata)
expressed_genes = which(apply(exprs_cqn, 1, mean) > 1)
exprs_cqn_filtered = exprs_cqn[expressed_genes, ]
heatmap.2(cor(exprs_cqn_filtered, method = "spearman"))

#Perform PCA analysis
pca_list = performPCA(exprs_cqn_filtered, design_matrix)
ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition_name, label = sample_id)) + 
  geom_text()

#Perform DEseq (This takes a really long time!)
dds_gene_meta = dplyr::select(gene_metadata,gene_id, gene_name, gene_biotype)
dds = DESeqDataSetFromMatrix(counts, design_filtered, ~SL1344 + IFNg + SL1344:IFNg)
dds = DESeq(dds)

#Extract IFNG results
ifng_results = results(dds, contrast = c("IFNg","primed","naive"))
ifng_table = ifng_results %>% data.frame() %>% dplyr::mutate(gene_id = rownames(ifng_results)) %>% 
  tbl_df() %>% 
  dplyr::left_join(dds_gene_meta, by = "gene_id") %>% 
  dplyr::arrange(padj)
ifng_up = dplyr::filter(ifng_table, padj < 0.01, log2FoldChange > 1, gene_biotype == "protein_coding") %>% arrange(-log2FoldChange)
ifng_down = dplyr::filter(ifng_table, padj < 0.01, log2FoldChange < -1, gene_biotype == "protein_coding") %>% arrange(log2FoldChange)
write.table(ifng_table, "results/SL1344/proteomics_comparison/DESeq2_IFNg_effect.txt", sep ="\t", quote = FALSE, row.names = FALSE)

#Extract Salmonella results
sl1344_results = results(dds, contrast = c("SL1344","infected","control"))
sl1344_table = sl1344_results %>% data.frame() %>% dplyr::mutate(gene_id = rownames(sl1344_results)) %>% 
  tbl_df() %>% 
  dplyr::left_join(dds_gene_meta, by = "gene_id") %>% 
  dplyr::arrange(padj)
sl1344_up = dplyr::filter(sl1344_table, padj < 0.01, log2FoldChange > 1, gene_biotype == "protein_coding") %>% arrange(-log2FoldChange)
sl1344_down = dplyr::filter(sl1344_table, padj < 0.01, log2FoldChange < -1, gene_biotype == "protein_coding") %>% arrange(log2FoldChange)
write.table(sl1344_table, "results/SL1344/proteomics_comparison/DESeq2_SL1344_effect.txt", sep ="\t", quote = FALSE, row.names = FALSE)

#Extract interaction results
interaction_results = results(dds, name = "SL1344infected.IFNgprimed")
interaction_table = interaction_results %>% data.frame() %>% dplyr::mutate(gene_id = rownames(interaction_results)) %>% 
  tbl_df() %>% 
  dplyr::left_join(dds_gene_meta, by = "gene_id") %>% 
  dplyr::arrange(padj)
interaction_up = dplyr::filter(interaction_table, padj < 0.01, log2FoldChange > 1, gene_biotype == "protein_coding") %>% arrange(-log2FoldChange)
interaction_down = dplyr::filter(interaction_table, padj < 0.01, log2FoldChange < -1, gene_biotype == "protein_coding") %>% arrange(log2FoldChange)
write.table(interaction_table, "results/SL1344/proteomics_comparison/DESeq2_IFNg_SL1344_interaction.txt", sep ="\t", quote = FALSE, row.names = FALSE)

#Export the counts data
design_filtered = design_filtered %>% dplyr::arrange(replicate, donor, condition)
write.table(design_filtered, "results/SL1344/proteomics_comparison/design_matrix.txt")
write.table(counts[,design_filtered$sample_id], "results/SL1344/proteomics_comparison/expression_counts.txt", sep = "\t", quote = FALSE)
write.table(exprs_cqn[,design_filtered$sample_id], "results/SL1344/proteomics_comparison/expression_quantile_normalised.txt", sep = "\t", quote = FALSE)

#Calculate mean expression per condition
exprs_cqn_mean = calculateMean(exprs_cqn, design_filtered, "condition")
exprs_cqn_mean = dplyr::mutate(exprs_cqn_mean, gene_id = rownames(exprs_cqn_mean)) %>% tbl_df()
pos_interaction = dplyr::semi_join(exprs_cqn_mean, interaction_up, by = "gene_id") %>% 
  dplyr::left_join(dds_gene_meta, by = "gene_id") %>%
  dplyr::mutate(salmonella_fc = C-A, interaction_fc = D-A)
synergestic_activation = dplyr::filter(pos_interaction, interaction_fc - salmonella_fc > 1, interaction_fc > 1)

#Plot examples of synergestic activation
design_filtered$condition_name = factor(design_filtered$condition_name, levels = c("naive","IFNg","SL1344", "IFNg+SL1344"))
plots = lapply(as.list(synergestic_activation$gene_id), plotGene, exprs_cqn, design_filtered, dds_gene_meta)
savePlots(plots, path = "results/SL1344/proteomics_comparison/",width = 5, height = 5)


#Compare to baseline
#Perform DEseq (This takes a really long time!)
dds_gene_meta = dplyr::select(gene_metadata,gene_id, gene_name, gene_biotype)
dds = DESeqDataSetFromMatrix(counts, design_filtered, ~condition_name)
dds = DESeq(dds)

#IFNG
ifng_results = results(dds, contrast = c("condition_name","IFNg","naive"))
ifng_table = ifng_results %>% data.frame() %>% dplyr::mutate(gene_id = rownames(ifng_results)) %>% 
  tbl_df() %>% 
  dplyr::left_join(dds_gene_meta, by = "gene_id") %>% 
  dplyr::arrange(padj)
write.table(ifng_table, "results/SL1344/proteomics_comparison/DESeq2_naive_vs_IFNg.txt", sep ="\t", quote = FALSE, row.names = FALSE)

#SL1344
sl1344_results = results(dds, contrast = c("condition_name","SL1344","naive"))
sl1344_table = sl1344_results %>% data.frame() %>% dplyr::mutate(gene_id = rownames(ifng_results)) %>% 
  tbl_df() %>% 
  dplyr::left_join(dds_gene_meta, by = "gene_id") %>% 
  dplyr::arrange(padj)
write.table(sl1344_table, "results/SL1344/proteomics_comparison/DESeq2_naive_vs_SL1344.txt", sep ="\t", quote = FALSE, row.names = FALSE)

#Both
both_results = results(dds, contrast = c("condition_name","IFNg.SL1344","naive"))
both_table = both_results %>% data.frame() %>% dplyr::mutate(gene_id = rownames(ifng_results)) %>% 
  tbl_df() %>% 
  dplyr::left_join(dds_gene_meta, by = "gene_id") %>% 
  dplyr::arrange(padj)
write.table(both_table, "results/SL1344/proteomics_comparison/DESeq2_naive_vs_IFNg+SL1344.txt", sep ="\t", quote = FALSE, row.names = FALSE)
