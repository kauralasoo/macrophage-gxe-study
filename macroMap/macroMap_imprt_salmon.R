library("readr")
library("tximport")
library("ggplot2")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping//")

sample_names = read.table("macrophage-gxe-study/macroMap/sample_names.txt", stringsAsFactors = FALSE)[,1]

#Construct design matrix
design = dplyr::data_frame(sample_id = sample_names) %>%
  tidyr::separate(sample_id, into = c("line_id", "treatment"), sep = "\\_FF\\_", remove = FALSE)

#Create a vector of file names
file_names = file.path("processed/macroMap/salmon/ensembl_87/", sample_names, "quant.sf")
names(file_names) = sample_names

#Import transcript metadata
transcript_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.transcript_data.rds")) %>%
  dplyr::rename(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id, gene_name = external_gene_name, chr = chromosome_name)
filtered_transcscript_data = readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.filtered.rds")

#Convert transcript data to suitable format for tximport
tx2gene = dplyr::select(transcript_data, gene_id, transcript_id, transcript_version) %>% 
  dplyr::mutate(TXNAME = paste(transcript_id, transcript_version, sep = ".")) %>% 
  dplyr::transmute(TXNAME, gene_id, transcript_id)

#Import gene-level abundances
gene_abundances = tximport(file_names, type = "salmon", tx2gene = tx2gene[,1:2], reader = read_tsv)
expressed_genes = gene_abundances$abundance[apply(gene_abundances$abundance, 1, mean) > 2,]

pca = performPCA(expressed_genes, design, n_pcs = 10)
ggplot(pca$pca_matrix, aes(x = PC1, y = PC2, label = treatment)) + geom_point() + geom_text()


rRNA_genes = dplyr::filter(transcript_data, gene_biotype == "rRNA")


