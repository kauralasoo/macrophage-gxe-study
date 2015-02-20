library("DESeq2")
library("dplyr")
library("devtools")
library("gplots")
library("ggplot2")
load_all("macrophage-gxe-study/seqUtils/")

#Load filtered gene metadata from disk
filtered_metadata = readRDS("annotations/biomart_transcripts.filtered.rds")
gene_data = dplyr::select(filtered_metadata, gene_id, gene_name, gene_biotype) %>% unique()

#Read count files from disk
sample_names = read.table("fastq/acLDL_samples.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,2]
data = read.table("results/acLDL/acLDL_combined_counts.txt", stringsAsFactors = FALSE, header = TRUE)

#Construct a design matrix
design = data.frame(sample_id = sample_names, treatment = c("ctrl", "acLDL", "ctrl", "acLDL", "ctrl", "acLDL"), stringsAsFactors = FALSE)
rownames(design) = design$sample_id

#Filter counts data
filtered_data = dplyr::filter(data, gene_id %in% gene_data$gene_id) #Remove pseudogenes
length_df = dplyr::select(filtered_data, gene_id, length) #Extract gene lengths
counts = dplyr::select(filtered_data, -gene_id, -length) #Extract counts matrix
rownames(counts) = filtered_data$gene_id

#Normalize data to transcripts per million (TPM)
tpm = calculateTPM(counts, lengths = length_df)
tpm_expressed = tpm[which(apply(tpm, 1, mean) > 2), ] #Keep only expressed genes

#Make heatmap of gene expression
cor_matrix = cor(log(tpm_expressed + 0.1, 2), method = "pearson")
heatmap.2(cor_matrix, margins = c(12,12))

#Perform PCA
tpm_z = zScoreNormalize(tpm_expressed)
pca_list = performPCA(tpm_z, design)
ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = treatment)) + geom_point()

#Peform differential expression analysis
dds = DESeqDataSetFromMatrix(count_matrix, design, ~treatment)
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
write.table(de_results, "results/acLDL/DE_genes.txt")
