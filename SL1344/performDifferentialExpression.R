library("DESeq2")
library("dplyr")
library("devtools")
library("tidyr")
load_all("macrophage-gxe-study/seqUtils/")

#Load data from disk
data = read.table("results//SL1344//combined_counts.txt", stringsAsFactors = FALSE)
filtered_metadata = readRDS("annotations/biomart_transcripts.filtered.rds")
gene_data = dplyr::select(filtered_metadata, gene_id, gene_name, gene_biotype) %>% unique()

#Filter genes
filtered_data = dplyr::filter(data, gene_id %in% gene_data$gene_id)
length_df = dplyr::select(filtered_data, gene_id, length)
counts = dplyr::select(filtered_data, -gene_id, -length)
rownames(counts) = filtered_data$gene_id

#TPM_normalize the data
tpm = calculateTPM(counts, lengths = length_df)

expressed_genes = which(apply(tpm, 1, mean) > 2)
tpm_expressed = tpm[expressed_genes, ]
heatmap.2(cor(log(tpm_expressed+0.1,2), method = "spearman"))

#Set up the design matrix
design = data.frame(sample_id = colnames(counts), stringsAsFactors = FALSE)
condition_names = data.frame(condition = c("A","B","C","D"), condition_name = c("naive", "IFNg", "SL1344", "IFNg+SL1344"), stringsAsFactors = FALSE)
design = design %>% 
  separate(sample_id, into = c("donor", "condition", "replicate"), extra = "drop", remove =FALSE) %>%
  mutate(replicate = ifelse(is.na(replicate), 1, replicate)) %>% #if no replicate in file name then set it to 1
  dplyr::mutate(replicate = ifelse(donor == "oarz", 1, replicate)) %>% #Oarz only has one replicate
  dplyr::left_join(condition_names, by = "condition") %>%
  dplyr::mutate(SL1344 = ifelse(condition %in% c("C","D"), "infected", "control")) %>%
  dplyr::mutate(IFNg = ifelse(condition %in% c("B","D"), "primed", "naive"))
rownames(design)= design$sample_id

dds = DESeqDataSetFromMatrix(counts, design, ~SL1344 + IFNg + SL1344:IFNg)
dds = DESeq(dds)


#Perform PCA analysis
tpm_z = zScoreNormalize(tpm_expressed)
pca = prcomp(t(tpm_z))
pca_res = as.data.frame(pca$x)
pca_matrix = as.data.frame(pca$x) %>% dplyr::mutate(sample_id = rownames(pca_res)) %>%
  dplyr::left_join(design, by = "sample_id")

