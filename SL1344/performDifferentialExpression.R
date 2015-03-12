library("DESeq2")
library("dplyr")
library("devtools")
library("tidyr")
library("ggplot2")
library("gplots")
library("cqn")
library("gProfileR")
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

#CQN normalize the data
exprs_cqn = calculateCQN(counts, gene_metadata)
expressed_genes = which(apply(exprs_cqn, 1, mean) > 1)
exprs_cqn_filtered = exprs_cqn[expressed_genes, ]
pdf("results/SL1344/diffExp/sample_heatmap.pdf", width = 10, height = 10)
heatmap.2(cor(exprs_cqn_filtered, method = "spearman"), trace = "none")
dev.off()

#Calculate mean expression in each condition
exprs_cqn_mean = calculateMean(exprs_cqn, design_matrix, "condition")
exprs_cqn_mean = dplyr::mutate(exprs_cqn_mean, gene_id = rownames(exprs_cqn_mean)) %>% tbl_df()

#Perform PCA analysis
pca_list = performPCA(exprs_cqn_filtered, design_matrix)
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
design_matrix$condition_name = factor(design_matrix$condition_name, levels = c("naive","IFNg","SL1344", "IFNg+SL1344"))
ifnb_plot = plotGene("ENSG00000171855",exprs_cqn, design_matrix, dds_gene_meta)
ggsave("results/SL1344/diffExp/IFNB1.pdf", plot = ifnb_plot, width = 5, height = 5)

