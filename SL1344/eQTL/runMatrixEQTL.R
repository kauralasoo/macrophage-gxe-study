library("MatrixEQTL")
library("ggplot2")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Load expresison dataset preapred previously by processExpressionData.R script
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/SL1344/array_genotypes.59_samples.vcfToMatrix.rds") #genotypes
gene_id_name_map = dplyr::select(expression_dataset$gene_metadata, gene_id, gene_name)

#Discard replicate samples
design = dplyr::filter(expression_dataset$design, !(donor == "fpdj")) %>% tbl_df() %>% #Remove all fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) %>% #Remove second ougl sample (dium)
  dplyr::filter(!(donor == "mijn")) #%>% #Remove mijn (wrong line from CGAP)
  #dplyr::filter(!(donor == "jorr")) #Very strong outlier in PEER analysis
sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))

#Set up expression data
exprs_cqn_all = expression_dataset$exprs_cqn[,design$sample_id]
expressed_genes = names(which(rowMeans(exprs_cqn_all) > 0)) #Set conservative threshold to expression level
exprs_cqn = exprs_cqn_all[expressed_genes,design$sample_id]

#Filter genotype data
sample_meta_filtered = dplyr::filter(sample_meta, condition == "A")
geno_data = extractSubset(sample_meta_filtered, vcf_file$genotypes, old_column_names = "genotype_id")

#Set up the genepos data.frame
genepos = dplyr::filter(expression_dataset$gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::transmute(geneid = gene_id, chr = chromosome_name, left = start_position, right = end_position) %>%
  as.data.frame()

#Set up expression data for each condition
condA_exp = extractSubset(dplyr::filter(sample_meta, condition == "A"), exprs_cqn)
condB_exp = extractSubset(dplyr::filter(sample_meta, condition == "B"), exprs_cqn)
condC_exp = extractSubset(dplyr::filter(sample_meta, condition == "C"), exprs_cqn)
condD_exp = extractSubset(dplyr::filter(sample_meta, condition == "D"), exprs_cqn)

#Construct covariate matrix
covariates = dplyr::mutate(sample_meta, sex = ifelse(gender == "male",1,0)) %>%
  dplyr::select(sample_id, sex, ng_ul_mean, diff_days)

#Import PEER factors
peer_factors = readRDS("results/SL1344/PEER/expressed_10_factors.rds")
peer_covariates = dplyr::left_join(covariates, peer_factors, by = "sample_id") %>% as.data.frame()
rownames(peer_covariates) = peer_covariates$sample_id
peer_covariates = t(peer_covariates[,-1])

#Set up covatiates for each conditon
condA_covs = extractSubset(dplyr::filter(sample_meta, condition == "A"), peer_covariates)
condB_covs = extractSubset(dplyr::filter(sample_meta, condition == "B"), peer_covariates)
condC_covs = extractSubset(dplyr::filter(sample_meta, condition == "C"), peer_covariates)
condD_covs = extractSubset(dplyr::filter(sample_meta, condition == "D"), peer_covariates)

#Run MatrixEQTL on each condition
condA_res_cov = runMatrixEQTL(condA_exp, geno_data, vcf_file$snpspos, genepos, covariates = condA_covs[1:7,])
condA_results = filterEQTLs(condA_res_cov$cis$eqtls, gene_id_name_map = gene_id_name_map, fdr_cutoff = 0.1)
dim(condA_results)
plot(condA_res_cov)

condB_res_cov = runMatrixEQTL(condB_exp, geno_data, vcf_file$snpspos, genepos, covariates = condB_covs[1:7,])
condB_results = filterEQTLs(condB_res_cov$cis$eqtls, gene_id_name_map = gene_id_name_map, fdr_cutoff = 0.1)
dim(condB_results)
plot(condB_res_cov)

condC_res_cov = runMatrixEQTL(condC_exp, geno_data, vcf_file$snpspos, genepos, covariates = condC_covs[1:7,])
condC_results = filterEQTLs(condC_res_cov$cis$eqtls, gene_id_name_map = gene_id_name_map, fdr_cutoff = 0.1)
dim(condC_results)
plot(condC_res_cov)

condD_res_cov = runMatrixEQTL(condD_exp, geno_data, vcf_file$snpspos, genepos, covariates = condD_covs[1:7,])
condD_results = filterEQTLs(condD_res_cov$cis$eqtls, gene_id_name_map = gene_id_name_map, fdr_cutoff = 0.1)
dim(condD_results)
plot(condD_res_cov)

#RUN matrixEQL
condA_res = runMatrixEQTL(condA_exp, geno_data, vcf_file$snpspos, genepos)
write.table(condA_res$cis$eqtls, "results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_A.txt", sep ="\t", quote = FALSE, row.names = FALSE)
condB_res = runMatrixEQTL(condB_exp, geno_data, vcf_file$snpspos, genepos)
write.table(condB_res$cis$eqtls, "results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_B.txt", sep ="\t", quote = FALSE, row.names = FALSE)
condC_res = runMatrixEQTL(condC_exp, geno_data, vcf_file$snpspos, genepos)
write.table(condC_res$cis$eqtls, "results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_C.txt", sep ="\t", quote = FALSE, row.names = FALSE)
condD_res = runMatrixEQTL(condD_exp, geno_data, vcf_file$snpspos, genepos)
write.table(condD_res$cis$eqtls, "results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_D.txt", sep ="\t", quote = FALSE, row.names = FALSE)

#### Save the data to disk for use with eqtlbma ####
#Expression data
write.table(condA_exp, "eqtlbma/input/exp_condition_A.txt", sep ="\t", quote = FALSE)
write.table(condB_exp, "eqtlbma/input/exp_condition_B.txt", sep ="\t", quote = FALSE)
write.table(condC_exp, "eqtlbma/input/exp_condition_C.txt", sep ="\t", quote = FALSE)
write.table(condD_exp, "eqtlbma/input/exp_condition_D.txt", sep ="\t", quote = FALSE)

#Gene positions
genepos_eqtlbma = dplyr::filter(expression_dataset$gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::transmute(chr = chromosome_name, left = start_position, right = end_position, geneid = gene_id, score = 1000, strand) %>%
  dplyr::mutate(strand = ifelse(strand == 1, "+","-")) %>%
  as.data.frame()
write.table(genepos_eqtlbma, "eqtlbma/input/gene_coords.bed", sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#SNP positions
snpspos_eqtlbma = dplyr::transmute(vcf_file$snpspos, chr, start  = pos-1, end = pos, snpid)
write.table(snpspos_eqtlbma, "eqtlbma/input/snp_coords.bed", sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Genotypes
write.table(geno_data, "eqtlbma/input/genotypes.txt", sep ="\t", quote = FALSE)


peer_factors = read.table("results/SL1344/PEER/expressed_B_10/factors.txt", sep =",")
peer_factors2 = read.table("results/SL1344/PEER/all_B_10/factors.txt", sep =",")

expression_list = list(naive = condA_exp, IFNg = condB_exp, SL1344 = condC_exp, IFNg_SL1344 = condC_exp)

