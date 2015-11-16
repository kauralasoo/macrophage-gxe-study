library("DESeq2")
library("dplyr")
library("devtools")
library("tidyr")
library("ggplot2")
library("gplots")
library("cqn")
library("gProfileR")
load_all("../seqUtils/")

#Load expression data from disk
#expression_list = readRDS("results/SL1344/combined_expression_data.rds")
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")

#Make heatmap
#Extract metadata for the heatmap command
meta = as.data.frame(dplyr::select(eqtl_data_list$sample_metadata, condition_name))
rownames(meta) = eqtl_data_list$sample_metadata$sample_id
#
pheatmap(cor(eqtl_data_list$exprs_cqn, method = "spearman"), annotation_row = meta, 
         show_rownames = FALSE, show_colnames = FALSE,
         filename = "results/SL1344/diffExp/sample_heatmap.pdf", width = 10, height = 8)

#Perform PCA analysis
pca_list = performPCA(eqtl_data_list$exprs_cqn, eqtl_data_list$sample_metadata)
pca_plot = ggplot(pca_list$pca_matrix, aes(x = PC1, y = PC2, color = condition_name, label = sample_id)) + 
  geom_text()
ggsave("results/SL1344/diffExp/PCA_of_gene_expression.pdf", plot = pca_plot, width = 8.5, height = 6.5)

#Perform DEseq (This takes a really long time!)
dds_gene_meta = dplyr::select(gene_metadata,gene_id, gene_name, gene_biotype)
dds = DESeqDataSetFromMatrix(counts, design_matrix, ~SL1344 + IFNg)
saveRDS(dds, "results/SL1344/DESeq2_results_full_data.rds")
dds = readRDS("results/SL1344/DESeq2_results_full_data.rds")

#IFNg results
ifng_results = results(dds, contrast = c("IFNg","primed","naive"))
ifng_results_filtered = filterDESeqResults(ifng_results, dds_gene_meta, biotype_filter = "protein_coding")

#Salmonella results
sl1344_results = results(dds, contrast = c("SL1344","infected","control"))
sl1344_results_filtered = filterDESeqResults(sl1344_results, dds_gene_meta, biotype_filter = "protein_coding")

#Interaction_genes
interaction_results = results(dds, name = "SL1344infected.IFNgprimed")
interaction_results_filtered = filterDESeqResults(interaction_results, dds_gene_meta, biotype_filter = "protein_coding")

#Iddentify syneregestically activated genes
pos_interaction = dplyr::semi_join(exprs_cqn_mean, interaction_results_filtered$up_genes, by = "gene_id") %>% 
  dplyr::left_join(dds_gene_meta, by = "gene_id") %>%
  dplyr::mutate(salmonella_fc = C-A, interaction_fc = D-A)
synergestic_activation = dplyr::filter(pos_interaction, interaction_fc - salmonella_fc > 1, interaction_fc > 1)


#Perform GO analysis
ifng_up_go = gprofiler(ifng_results_filtered$up_genes$gene_id, organism = "hsapiens") %>% tbl_df() %>% arrange(p.value)
sl1344_up_go = gprofiler(sl1344_results_filtered$up_genes$gene_id, organism = "hsapiens") %>% tbl_df() %>% arrange(p.value)
dplyr::filter(sl1344_up_go, term.size < 300) %>% dplyr::select(p.value, term.name, term.size) %>% View()

go_results = list(ifng_up_go = ifng_up_go, sl1344_up_go = sl1344_up_go)
saveRDS(go_results, "results/SL1344/diffExp/GO_results.rds")

#Make plot of IFNb
design_matrix = expression_list$design
design_matrix$condition_name = factor(design_matrix$condition_name, levels = c("naive","IFNg","SL1344", "IFNg_SL1344"))
ifnb_plot = plotGene("ENSG00000171855",expression_list$exprs_cqn, design_matrix, expression_list$gene_metadata)
ggsave("results/SL1344/diffExp/IFNB1.pdf", plot = ifnb_plot, width = 5, height = 5)

#Make some plots for subho
subho_genes = read.table("../LPS/subho_genes.txt", stringsAsFactors = FALSE)[,1]
plots = lapply(as.list(subho_genes), plotGene, expression_list$exprs_cqn, design_matrix, expression_list$gene_metadata)
savePlots(plots, "subho_plots/", 6, 6)
