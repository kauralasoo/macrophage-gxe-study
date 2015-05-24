library("devtools")
library("cqn")
load_all("macrophage-gxe-study/seqUtils/")

#Import raw read counts 
data = read.table("results/SL1344/SL1344_basic_counts.txt", stringsAsFactors = FALSE, header = TRUE)
#There seems to be a sample swap between coio_C and coio_D, fix that
indexes = c(which(colnames(data) == "coio_C"), which(colnames(data) == "coio_D"))
colnames(data)[indexes] = c("coio_D", "coio_C")
#Also, all vorx and zuta samples seem to have been swapped. Fix it.
samples = colnames(data)
vorx_samples = c(which(samples == "vorx_A"),
                 which(samples == "vorx_B"),
                 which(samples == "vorx_C"),
                 which(samples == "vorx_D"))
zuta_samples = c(which(samples == "zuta_A"),
                 which(samples == "zuta_B"),
                 which(samples == "zuta_C"),
                 which(samples == "zuta_D"))
samples[vorx_samples] = c("zuta_A","zuta_B","zuta_C","zuta_D")
samples[zuta_samples] = c("vorx_A","vorx_B","vorx_C","vorx_D")
colnames(data) = samples

#Load transcript metadata
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")
metadata = readRDS("annotations/Homo_sapiens.GRCh38.79.transcript_data.rds")

#Filter annotations
filtered_metadata = dplyr::filter(metadata,transcript_gencode_basic == "GENCODE basic") %>% 
  dplyr::select(ensembl_gene_id, gene_biotype, chromosome_name, percentage_gc_content,external_gene_name) %>% 
  unique() %>%
  dplyr::filter(gene_biotype %in% valid_gene_biotypes, chromosome_name %in% valid_chromosomes) %>%
  dplyr::rename(gene_id = ensembl_gene_id, gene_name = external_gene_name)

#Add gene length to the metadata
length_df = dplyr::select(data, gene_id, length)
gene_metadata = dplyr::left_join(filtered_metadata, length_df, by = "gene_id") %>% tbl_df()

#Process counts
filtered_data = dplyr::filter(data, gene_id %in% gene_metadata$gene_id)
counts = dplyr::select(filtered_data, -gene_id, -length)
rownames(counts) = filtered_data$gene_id
counts = counts[gene_metadata$gene_id,] #Reoder counts

#Construct a design matrix from the sample names
design_matrix = constructDesignMatrix_SL1344(sample_ids = colnames(counts))

#cqn_normalize_data
exprs_cqn = calculateCQN(counts, gene_metadata)

results_list = list(
  exprs_counts = counts,
  exprs_cqn = exprs_cqn,
  design = design_matrix,
  gene_metadata = gene_metadata)

saveRDS(results_list, "results/SL1344/combined_expression_data.rds")
