library("MatrixEQTL")
library("ggplot2")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Import eqtl data set from disk
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Use the first 7 covariates
eqtl_data_list$covariates_list = lapply(eqtl_data_list$covariates_list, function(x){x[1:7,]})

#Reorder columns in genpos to fit matrixEQTL requirements
eqtl_data_list$genepos = dplyr::select(eqtl_data_list$genepos, geneid, chr, left, right)

#Run MatrixEQTL on each condition
condA_res_cov = runMatrixEQTL(eqtl_data_list$exprs_cqn_list$naive, eqtl_data_list$genotypes, eqtl_data_list$snpspos, 
                              eqtl_data_list$genepos, covariates = eqtl_data_list$covariates_list$naive)
condA_results = filterEQTLs(condA_res_cov$cis$eqtls, gene_id_name_map = gene_id_name_map, fdr_cutoff = 0.1)
dim(condA_results)
plot(condA_res_cov)

condB_res_cov = runMatrixEQTL(eqtl_data_list$exprs_cqn_list$IFNg, eqtl_data_list$genotypes, eqtl_data_list$snpspos, 
                              eqtl_data_list$genepos, covariates = eqtl_data_list$covariates_list$IFNg)
condB_results = filterEQTLs(condB_res_cov$cis$eqtls, gene_id_name_map = gene_id_name_map, fdr_cutoff = 0.1)
dim(condB_results)
plot(condB_res_cov)

condC_res_cov = runMatrixEQTL(eqtl_data_list$exprs_cqn_list$SL1344, eqtl_data_list$genotypes, eqtl_data_list$snpspos, 
                              eqtl_data_list$genepos, covariates = eqtl_data_list$covariates_list$SL1344)
condC_results = filterEQTLs(condC_res_cov$cis$eqtls, gene_id_name_map = gene_id_name_map, fdr_cutoff = 0.1)
dim(condC_results)
plot(condC_res_cov)

condD_res_cov = runMatrixEQTL(eqtl_data_list$exprs_cqn_list$IFNg_SL1344, eqtl_data_list$genotypes, eqtl_data_list$snpspos, 
                              eqtl_data_list$genepos, covariates = eqtl_data_list$covariates_list$IFNg_SL1344)
condD_results = filterEQTLs(condD_res_cov$cis$eqtls, gene_id_name_map = gene_id_name_map, fdr_cutoff = 0.1)
dim(condD_results)
plot(condD_res_cov)
