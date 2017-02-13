library("tximport")
library("readr")
load_all("macrophage-gxe-study/housekeeping//")

#Import sample names
sample_names = read.table("macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_all.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

#Import transcript metadata
transcript_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.transcript_data.rds"))

#Create a vector of file names
file_names = file.path("processed/acLDL/salmon/ensembl_87/", sample_names, "quant.sf")
names(file_names) = sample_names

#Convert transcript data to suitable format for tximport
tx2gene = dplyr::select(transcript_data, ensembl_gene_id, ensembl_transcript_id, transcript_version) %>% 
  dplyr::mutate(TXNAME = paste(ensembl_transcript_id, transcript_version, sep = ".")) %>% 
  dplyr::transmute(TXNAME,GENEID = ensembl_gene_id)

#Import gene-level abundances
txi = tximport(file_names, type = "salmon", tx2gene = tx2gene, reader = read_tsv)

#Perform PCA
design = constructDesignMatrix_acLDL(sample_names)
filtered_genes = txi$abundance[rowMeans(txi$abundance) > 2,]
log2_tpm = log(filtered_genes + 0.01, 2)
pca = performPCA(log2_tpm, design, n_pcs = 5)
ggplot(a$pca_matrix, aes(x = PC1, y = PC2, color = condition_name, label = sample_id)) + geom_text()


ggplot(dplyr::filter(pca$pca_matrix, donor %in% c("hehd", "fikt", "qolg", "coyi","nukw")), aes(x = PC1, y = PC2, color = condition_name, label = sample_id)) + geom_text()
