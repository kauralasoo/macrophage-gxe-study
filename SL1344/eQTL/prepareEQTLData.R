library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("SNPRelate")


#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")

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


### FastQTL ####
#Export expression data for fastQTL analysis
eqtl_dataset = readRDS("results/SL1344/eqtl_data_list.rds")

#Extract data from list
donor_genotype_map = dplyr::filter(eqtl_dataset$sample_metadata, condition == "A") %>% dplyr::select(donor, genotype_id)
genepos = dplyr::select(eqtl_dataset$genepos, chr, left, right, geneid)
colnames(genepos)[1] = "#chr"

#Make an updated list of covariates
fastqtl_covariates_list = lapply(eqtl_dataset$covariates_list, function(cov_mat, donor_genotype_map){
  res = extractSubset(donor_genotype_map, as.data.frame(cov_mat), "donor", "genotype_id") %>% 
    dplyr::mutate(id = rownames(cov_mat)) %>% 
    dplyr::select(id, everything())
  return(res)
}, donor_genotype_map) %>%
  lapply(., function(x){x[1:7,]}) #Keep the first 7 covariates

#Make an updated list of gene expression
fastqtl_expression_list = lapply(eqtl_dataset$exprs_cqn_list, function(exp_mat, donor_genotype_map, genepos){
  res = extractSubset(donor_genotype_map, as.data.frame(exp_mat), "donor","genotype_id") %>%
    dplyr::mutate(geneid = rownames(exp_mat)) %>%
    dplyr::select(geneid, everything()) %>%
    dplyr::left_join(genepos, ., by = "geneid")
  return(res)
}, donor_genotype_map, genepos)

#Save matrices to disk
saveFastqtlMatrices(fastqtl_covariates_list, "results/SL1344/fastqtl/input/", file_suffix = "covariates")
saveFastqtlMatrices(fastqtl_expression_list, "results/SL1344/fastqtl/input/", file_suffix = "expression")
write.table(donor_genotype_map$genotype_id, "results/SL1344/fastqtl/input/genotype_list.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:200), n = 200)
write.table(chunks_matrix, "results/ATAC/fastqtl/input/chunk_table.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")

#### eigenMT ####
SNPRelate::snpgdsVCF2GDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.vcf.gz", 
                         "genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.gds", method = "copy.num.of.ref")
vcf_file = gdsToMatrix("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.gds")
saveRDS(vcf_file, "genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

eigenMTExportGenotypes(vcf_file, combined_expression_data$gene_metadata, "results/SL1344/eigenMT/")

#Save SNP positions
snp_pos_df = vcf_file$snpspos %>% dplyr::rename(snp = snpid, chr_snp = chr)
write.table(snp_pos_df, "results/SL1344/eigenMT/snp_positions.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Save genotypes
genotypes = dplyr::mutate(as.data.frame(vcf_file$genotypes), ID = rownames(vcf_file$genotypes)) %>% dplyr::select(ID, everything())
write.table(genotypes, "results/SL1344/eigenMT/genotypes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Save gene positions
gene_data = dplyr::transmute(combined_expression_data$gene_metadata, gene_id, chrom_probe = chromosome_name, s1 = start_position, s2 = end_position) %>%
  dplyr::arrange(chrom_probe, s1)
write.table(gene_data, "results/SL1344/eigenMT/gene_positions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
