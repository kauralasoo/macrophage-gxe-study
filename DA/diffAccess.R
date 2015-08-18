library("modules")
library("dplyr")
library("DESeq2")
unload(hk)
hk = import("macrophage-chromatin/housekeeping/constructDesignMatrices")

sample_names = read.table("macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]
design_matrix = hk$constructDesignMatrix_ATAC(sample_names)

#Test 2-way model with DEseq2
diff_design = dplyr::filter(design_matrix, donor %in% c("eofe","cicb","vass","bima")) %>% as.data.frame()
rownames(diff_design) = diff_design$sample_id

#Import counts
counts = readRDS("results/ATAC//ATAC_combined_counts.rds")
rownames(counts) = counts$gene_id
diff_counts = counts[,rownames(diff_design)]

dds = DESeqDataSetFromMatrix(diff_counts, diff_design, ~SL1344 + IFNg + IFNg:SL1344)
dds = DESeq(dds)

#Extract Salmonella results
sl1344_results = results(dds, name = "SL1344_infected_vs_control")
sl1344_results$peak_id = rownames(sl1344_results)
sl1344_results_filtered = dplyr::filter(as.data.frame(sl1344_results), padj < 0.05, abs(log2FoldChange) > 2) %>% 
  dplyr::arrange(padj)

#Extract IFNg results
ifng_results = results(dds, name = "IFNg_primed_vs_naive")
ifng_results$peak_id = rownames(ifng_results)
ifng_results_filtered = dplyr::filter(as.data.frame(ifng_results), padj < 0.05, abs(log2FoldChange) > 2) %>% 
  dplyr::arrange(padj)

#Extract IFNg results
interaction_results = results(dds, name = "SL1344infected.IFNgprimed")
interaction_results$peak_id = rownames(interaction_results)
interaction_results_filtered = dplyr::filter(as.data.frame(interaction_results), padj < 0.05, abs(log2FoldChange) > 2) %>% 
  dplyr::arrange(padj)