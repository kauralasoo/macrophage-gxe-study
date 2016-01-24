library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")


#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$gene_metadata = dplyr::rename(combined_expression_data$gene_metadata, 
                                                       chr = chromosome_name, start = start_position, end = end_position)

#Extract separate lists for each condition
condition_names = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
rna_conditions = lapply(condition_names, extractConditionFromExpressionList, combined_expression_data)

#Rename column names to genotype ids
rna_conditions_renamed = lapply(rna_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#### RASQUAL ####
rasqual_input_folder = "results/SL1344/rasqual/input/"
exportDataForRasqual(rna_conditions_renamed, rasqual_input_folder)

#Import exon coordinates from disk
union_exon_coords = read.table("annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.filtered_genes.txt", 
                               stringsAsFactors = FALSE, header = TRUE) %>% dplyr::rename(chr = chromosome_name)
filtered_coords = dplyr::semi_join(union_exon_coords, rna_conditions_renamed$naive$gene_metadata, by = "gene_id")
snp_coords = readr::read_delim("genotypes/SL1344/imputed_20151005/imputed.86_samples.snp_coords.txt", 
                               delim = "\t", col_types = "cdc", col_names = c("chr","pos","snp_id"))
#snp_coords = readr::read_delim("genotypes/SL1344/imputed_20151005/imputed.86_samples.snp_coords.INFO_07.txt", 
#                               delim = "\t", col_types = "cdc", col_names = c("chr","pos","snp_id"))

#Count the numer of overlapping SNPs
exon_df = countSnpsOverlapingExons(filtered_coords, snp_coords, cis_window = 500000) %>% dplyr::arrange(chromosome_name, range_start)
write.table(exon_df, file.path(rasqual_input_folder, "gene_snp_count_500kb.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

exon_df = countSnpsOverlapingExons(filtered_coords, snp_coords, cis_window = 100000) %>% dplyr::arrange(chromosome_name, range_start)
write.table(exon_df, file.path(rasqual_input_folder, "gene_snp_count_100kb.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

#Construct batches
chr11_batches = dplyr::filter(exon_df, chromosome_name == "11") %>%
  seqUtils::rasqualConstructGeneBatches(10)
write.table(chr11_batches, file.path(rasqual_input_folder, "chr11_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#Split genes into batches based on how many cis and feature SNPs they have
batches = rasqualOptimisedGeneBatches(exon_df, c(20,8,3,1))
write.table(batches, file.path(rasqual_input_folder, "gene_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

### FastQTL ####
fastqtl_input_folder = "results/SL1344/fastqtl/input/"
exportDataForFastQTL(rna_conditions_renamed, fastqtl_input_folder)

#Export chr11 data only
chr11_genes = dplyr::filter(combined_expression_data$gene_metadata, chr == 11)$gene_id
chr11_data = lapply(rna_conditions_renamed, extractGenesFromExpressionList, chr11_genes)
chr11_input_folder = "results/SL1344/fastqtl/input_chr11/"
exportDataForFastQTL(chr11_data, chr11_input_folder)

#### eigenMT ####
#Export genotype data
chromosome_list = scan("macrophage-gxe-study/data/sample_lists/chromosome_list.txt", what = "char")
eigenMTExportGenotypesByChr(chromosome_list, "genotypes/SL1344/imputed_20151005/chromosomes_INFO_07//",
                            "results/SL1344/eigenMT/input/", "chr_")
#Export gene metadata
eigenMTExportGeneMetadata(combined_expression_data$gene_metadata, "results/SL1344/eigenMT/input/")
