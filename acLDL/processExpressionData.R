library("devtools")
library("cqn")
library("dplyr")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping//")

#Import raw read counts 
data = read.table("results/acLDL/acLDL_basic_counts.txt", stringsAsFactors = FALSE, header = TRUE)

#Import metadata
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")
metadata = readRDS("annotations/Homo_sapiens.GRCh38.79.transcript_data.rds")

#Filter annotations
filtered_metadata = dplyr::filter(metadata,transcript_gencode_basic == "GENCODE basic") %>% 
  dplyr::select(ensembl_gene_id, gene_biotype, chromosome_name, strand, start_position, end_position, 
                percentage_gc_content, external_gene_name) %>% 
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
design_matrix = constructDesignMatrix_acLDL(sample_ids = colnames(counts))

#cqn_normalize_data
exprs_cqn = calculateCQN(counts, gene_metadata)

results_list = list(
  exprs_counts = counts,
  exprs_cqn = exprs_cqn,
  design = design_matrix,
  gene_metadata = gene_metadata)

saveRDS(results_list, "results/acLDL/acLDL_combined_expression_data.rds")

