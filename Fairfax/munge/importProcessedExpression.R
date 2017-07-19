library("dplyr")
library("readr")
library("biomaRt")
library("rtracklayer")
library("SummarizedExperiment")
library("devtools")
load_all("../seqUtils/")

#Import all matrices
cd14 = readr::read_tsv("../../datasets/Fairfax_2014/processed_expression/CD14.47231.414.b.txt.gz", col_names = TRUE)
ifn = readr::read_tsv("../../datasets/Fairfax_2014/processed_expression/IFN.47231.367.b.txt.gz", col_names = TRUE)
lps2 = readr::read_tsv("../../datasets/Fairfax_2014/processed_expression/LPS2.47231.261.b.txt.gz", col_names = TRUE)
lps24 = readr::read_tsv("../../datasets/Fairfax_2014/processed_expression/LPS24.47231.322.b.txt.gz", col_names = TRUE)

#Fetch gene ids from Ensembl
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host='www.ensembl.org')
gene_metadata = getBM(attributes=c("illumina_humanht_12_v4", "ensembl_gene_id", "chromosome_name", "start_position",
                                   "end_position", "strand","gene_biotype", "external_gene_name"), mart = ensembl)

valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")
valid_chomosome_names = c(as.character(c(1:22)), "X","Y")

#Filter to uniquely mapping probes
filtered_metadata = tbl_df(gene_metadata) %>% 
  dplyr::filter(illumina_humanht_12_v4 != "") %>% 
  dplyr::arrange(illumina_humanht_12_v4) %>% 
  dplyr::rename(probe_id = illumina_humanht_12_v4, 
                gene_id = ensembl_gene_id,
                chr = chromosome_name,
                gene_start = start_position,
                gene_end = end_position,
                gene_name = external_gene_name) %>%
  dplyr::group_by(probe_id) %>%
  dplyr::mutate(gene_count = length(gene_id)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(gene_count == 1) %>%
  dplyr::select(-gene_count) %>%
  dplyr::filter(gene_biotype %in% valid_gene_biotypes)

#Import Ensembl_74 gff file with GRCh37 coordinates
ensembl_74_gff = readr::read_tsv("../../annotations/GRCh37/Ensembl_74/Homo_sapiens.GRCh37.74.gff3", skip = 1, col_names = FALSE)
ensembl74_meta = dplyr::filter(ensembl_74_gff, X3 == "gene") %>% 
  tidyr::separate(X9, c("F","S","gene_id"), sep = ";Name=") %>% 
  dplyr::transmute(gene_id, chr = X1, gene_start = X4, gene_end = X5, strand = X7)

#Add new gene coordinates to metadata
metadata_new_coords = dplyr::select(filtered_metadata, -chr, -gene_start, -gene_end, -strand) %>% 
  dplyr::left_join(ensembl74_meta, by = "gene_id") %>% dplyr::filter(!is.na(gene_start)) %>%
  dplyr::filter(chr %in% valid_chomosome_names)

final_gene_metadata = as.data.frame(metadata_new_coords)
rownames(final_gene_metadata) = final_gene_metadata$probe_id

#Make a single expression matrix
cd14_mat = as.matrix(cd14[,-1])
rownames(cd14_mat) = cd14$PROBE_ID
colnames(cd14_mat) = paste("CD14", colnames(cd14_mat), sep = "_")

lps2_mat = as.matrix(lps2[,-1])
rownames(lps2_mat) = lps2$PROBE_ID
colnames(lps2_mat) = paste("LPS2", colnames(lps2_mat), sep = "_")

lps24_mat = as.matrix(lps24[,-1])
rownames(lps24_mat) = lps24$PROBE_ID
colnames(lps24_mat) = paste("LPS24", colnames(lps24_mat), sep = "_")

ifn_mat = as.matrix(ifn[,-1])
rownames(ifn_mat) = ifn$PROBE_ID
colnames(ifn_mat) = paste("IFN", colnames(ifn_mat), sep = "_")

matrix_list = list(cd14_mat[final_gene_metadata$probe_id,], lps2_mat[final_gene_metadata$probe_id,], 
                   lps24_mat[final_gene_metadata$probe_id,], ifn_mat[final_gene_metadata$probe_id,])
expression_matrix = purrr::reduce(matrix_list, cbind)

#Construct sample metadata
sample_metadata = data_frame(sample_id = colnames(expression_matrix)) %>% 
  tidyr::separate(sample_id, c("condition_name", "donor_id"), sep = "_", remove = FALSE) %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(condition_count = length(condition_name)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(present_in_all = ifelse(condition_count == 4, TRUE, FALSE))

#Count IFN diff samples
ifn_present = dplyr::filter(sample_metadata, condition_name %in% c("CD14", "IFN")) %>% 
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(condition_count = length(condition_name)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(condition_count == 2) %>%
  dplyr::mutate(present_in_IFN = TRUE) %>%
  dplyr::select(sample_id, present_in_IFN)

final_sample_metadata = dplyr::left_join(sample_metadata, ifn_present, by = "sample_id") %>%
  dplyr::mutate(present_in_IFN = ifelse(is.na(present_in_IFN), FALSE, TRUE)) %>%
  dplyr::mutate(genotype_id = paste(donor_id, donor_id, sep = "_")) %>%
  as.data.frame()
rownames(final_sample_metadata) = final_sample_metadata$sample_id

#Create summarised experiment matrix
se = SummarizedExperiment::SummarizedExperiment(
  assays = list(exprs = expression_matrix), 
  colData = final_sample_metadata, 
  rowData = final_gene_metadata)
saveRDS(se, "results/Fairfax/expression_data.SummarizedExperiment.rds")

