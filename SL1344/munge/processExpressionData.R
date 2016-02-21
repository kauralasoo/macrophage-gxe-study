library("devtools")
library("cqn")
library("dplyr")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Import raw read counts 
data = read.table("results/SL1344/SL1344_basic_counts.txt", stringsAsFactors = FALSE, header = TRUE)

#Load transcript metadata
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")
metadata = readRDS("annotations/Homo_sapiens.GRCh38.79.transcript_data.rds")
sample_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_salmonella_metadata.rds")

#Filter annotations
filtered_metadata = dplyr::filter(metadata,transcript_gencode_basic == "GENCODE basic") %>% 
  dplyr::select(ensembl_gene_id, gene_biotype, chromosome_name, strand, start_position, end_position, 
                percentage_gc_content, external_gene_name) %>% 
  unique() %>%
  dplyr::filter(gene_biotype %in% valid_gene_biotypes, chromosome_name %in% valid_chromosomes) %>%
  dplyr::rename(gene_id = ensembl_gene_id, gene_name = external_gene_name)

#Add exon start and end coordinates
union_exon_coords = read.table("annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.filtered_genes.txt", 
                               stringsAsFactors = FALSE, header = TRUE) %>% dplyr::rename(chr = chromosome_name)
filtered_coords = dplyr::semi_join(union_exon_coords, filtered_metadata, by = "gene_id") %>%
  dplyr::select(gene_id, exon_starts, exon_ends)

#Add gene length to the metadata
length_df = dplyr::select(data, gene_id, length)
gene_metadata = dplyr::left_join(filtered_metadata, length_df, by = "gene_id") %>% tbl_df() %>%
  dplyr::rename(chr = chromosome_name, start = start_position, end = end_position) %>%
  dplyr::left_join(filtered_coords, by = "gene_id")

#Process counts
filtered_data = dplyr::filter(data, gene_id %in% gene_metadata$gene_id)
counts = dplyr::select(filtered_data, -gene_id, -length)
rownames(counts) = filtered_data$gene_id
counts = counts[gene_metadata$gene_id,] #Reoder counts

#Filter out some samples and discard replicates from the design_matrix
design_matrix = constructDesignMatrix_SL1344(sample_ids = colnames(counts)) %>% #Construct a design matrix from the sample names
  dplyr::filter(!(donor == "fpdj")) %>% tbl_df() %>% #Remove both fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) %>% #Remove second ougl sample (dium)
  dplyr::filter(!(donor == "mijn")) %>% #Remove mijn (wrong line from CGAP) 
  dplyr::arrange(donor)

#Add metadata to the design matrix
sample_meta = dplyr::left_join(design_matrix, sample_metadata, by = c("donor","replicate"))

#Filter counts to include only samples that are in the design matrix
counts = counts[, design_matrix$sample_id]

#cqn_normalize_data
exprs_cqn = calculateCQN(counts, gene_metadata)

#Normalize data using TPM
exprs_tpm = calculateTPM(as.matrix(counts), as.data.frame(gene_metadata), fragment_length = 250)

#Calculate normalization factors
norm_factors = calculateNormFactors(counts, method = "RLE")

#Combine everything into a list
results_list = list(
  counts = counts,
  cqn = exprs_cqn,
  tpm = exprs_tpm,
  norm_factors = norm_factors,
  sample_metadata = sample_meta,
  gene_metadata = gene_metadata)

saveRDS(results_list, "results/SL1344/combined_expression_data.rds")
