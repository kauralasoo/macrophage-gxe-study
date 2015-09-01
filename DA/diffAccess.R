library("modules")
library("dplyr")
library("DESeq2")
source("macrophage-chromatin/housekeeping/constructDesignMatrices.r")

#Construct design matrix
sample_names = read.table("macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
sample_names = sample_names[grepl("H3K27Ac",sample_names)]
design = constructDesignMatrix_Ivashkiv(sample_names)
rownames(design) = design$sample_id

#Import peak counts
counts = readRDS("results/Ivashkiv/H3K27Ac_combined_counts.rds")
rownames(counts) = counts$gene_id
diff_counts = counts[,rownames(design)]

#Run DEseq
dds = DESeqDataSetFromMatrix(diff_counts, design, ~condition)
dds = DESeq(dds)

ifng_results = results(dds, contrast = c("condition","naive","IFNg"))
ifng_results$peak_id = rownames(ifng_results)
ifng_results_filtered = dplyr::filter(as.data.frame(ifng_results), padj < 0.05, abs(log2FoldChange) > 1) %>% 
  dplyr::arrange(padj)


lps_results = results(dds, contrast = c("condition","naive","LPS"))
lps_results$peak_id = rownames(lps_results)
lps_results_filtered = dplyr::filter(as.data.frame(lps_results), padj < 0.05, abs(log2FoldChange) > 1) %>% 
  dplyr::arrange(padj)

both_results = results(dds, contrast = c("condition","naive","IFNg_LPS"))
both_results$peak_id = rownames(both_results)
both_results_filtered = dplyr::filter(as.data.frame(both_results), padj < 0.05, abs(log2FoldChange) > 1) %>% 
  dplyr::arrange(padj)