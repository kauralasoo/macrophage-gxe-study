library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data.rds")

#Filter genes by expression level
mean_expression = calculateMean(combined_expression_data$cqn, as.data.frame(combined_expression_data$sample_metadata), "condition_name")
expressed_genes = names(which(apply(mean_expression, 1, max) > 0))
not_Y_genes = dplyr::filter(combined_expression_data$gene_metadata, !(chromosome_name %in% c("MT","Y")))$gene_id
keep_genes = intersect(expressed_genes, not_Y_genes)
rna_expressed = extractGenesFromExpressionList(combined_expression_data, keep_genes)

#Extract separate lists for each condition
condition_names = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
rna_conditions = lapply(condition_names, extractConditionFromExpressionList, rna_expressed)

#Rename column names to genotype ids
rna_conditions_renamed = lapply(rna_conditions, renameMatrixColumnsInExpressionList, "sample_id", "genotype_id")

#### export data for RASQUAL ####
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

#Construct batches
chr11_batches = dplyr::filter(exon_df, chromosome_name == "11") %>%
  seqUtils::rasqualConstructGeneBatches(10)
write.table(chr11_batches, file.path(rasqual_input_folder, "chr11_batches.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#### COVARIATES ####
#Construct covariate matrix
covariates = dplyr::mutate(sample_meta, sex = ifelse(gender == "male",1,0)) %>%
  dplyr::select(sample_id, sex, ng_ul_mean, diff_days)

#Save expression data for PEER and run PEER outside of R
savePEERData(exprs_cqn_list, "results/SL1344/PEER/input/")

#Run PEER .....

#Import PEER results back into R
cond_A_factors = importPEERFactors("results/SL1344/PEER/naive_10/factors.txt", cond_A_design)
cond_B_factors = importPEERFactors("results/SL1344/PEER/IFNg_10/factors.txt", cond_B_design)
cond_C_factors = importPEERFactors("results/SL1344/PEER/SL1344_10/factors.txt", cond_C_design)
cond_D_factors = importPEERFactors("results/SL1344/PEER/IFNg_SL1344_10/factors.txt", cond_D_design)
peer_factors = rbind(cond_A_factors,cond_B_factors, cond_C_factors, cond_D_factors) %>%
  dplyr::semi_join(design, by = "sample_id")

#Merge PEER facotrs with covariates
peer_covariates = dplyr::left_join(covariates, peer_factors, by = "sample_id") %>% as.data.frame()
rownames(peer_covariates) = peer_covariates$sample_id
peer_covariates = t(peer_covariates[,-1])

#Set up covatiates for each conditon
condA_covs = extractSubset(dplyr::filter(sample_meta, condition == "A"), peer_covariates)
condB_covs = extractSubset(dplyr::filter(sample_meta, condition == "B"), peer_covariates)
condC_covs = extractSubset(dplyr::filter(sample_meta, condition == "C"), peer_covariates)
condD_covs = extractSubset(dplyr::filter(sample_meta, condition == "D"), peer_covariates)
covariates_list = list(naive = condA_covs, IFNg = condB_covs, SL1344 = condC_covs, IFNg_SL1344 = condD_covs)

#### Dataset object ####
#Construct a single dataset object that can be reused by many different analysis
eqtl_dataset = list(
  exprs_cqn_list = exprs_cqn_list,
  covariates_list = covariates_list,
  genepos = genepos,
  gene_metadata = expression_dataset$gene_metadata,
  sample_metadata = sample_meta,
  covariates = t(peer_covariates),
  exprs_cqn = exprs_cqn
  )
saveRDS(eqtl_dataset, "results/SL1344/eqtl_data_list.rds")

#### eqtlbma ####
#Export eqtl dataset for analysis with eqtlbma
#Use only the first 7 covariates
eqtlbma_dataset = eqtl_dataset
eqtlbma_dataset$covariates_list = lapply(eqtl_dataset$covariates_list, function(x){x[1:7,]})
#Save data for eqtlbma
saveEqtlbmaData(eqtlbma_dataset, "results/SL1344/eqtlbma/input/", 
                project_root = "/nfs/users/nfs_k/ka8/group-scratch/kaur/projects/macrophage-gxe-study/")

#### RASQUAL ####
#export eqtl dataset for analysis with rasqual


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

