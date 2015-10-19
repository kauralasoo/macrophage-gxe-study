library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Import imputed genotype data from disk
#SNPRelate::snpgdsVCF2GDS("genotypes/SL1344/imputed_20151005/imputed.59_samples.snps_indels.INFO_08.vcf.gz", "genotypes/SL1344/imputed_20151005/imputed.59_samples.snps_indels.INFO_08.gds",method = "copy.num.of.ref")
vcf_file = gdsToMatrix("genotypes/SL1344/imputed_20151005/imputed.59_samples.snps_indels.INFO_08.gds")

#Import genotype data from the VCF file
#vcf_file = vcfToMatrix("genotypes/SL1344/array_genotypes.59_samples.imputed.uniq.vcf", "GRCh38")
#saveRDS(vcf_file, "genotypes/SL1344/array_genotypes.59_samples.imputed.vcfToMatrix.rds")

#Load the raw eQTL dataset
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
gene_id_name_map = dplyr::select(expression_dataset$gene_metadata, gene_id, gene_name)

#Filter out some samples and discard replicates
design = dplyr::filter(expression_dataset$design, !(donor == "fpdj")) %>% tbl_df() %>% #Remove both fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) %>% #Remove second ougl sample (dium)
  dplyr::filter(!(donor == "mijn")) #%>% #Remove mijn (wrong line from CGAP)
sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))

#Construct separate design matrices
cond_A_design = dplyr::filter(sample_meta, condition == "A")
cond_B_design = dplyr::filter(sample_meta, condition == "B")
cond_C_design = dplyr::filter(sample_meta, condition == "C")
cond_D_design = dplyr::filter(sample_meta, condition == "D")
design_list = list(naive = cond_A_design, IFNg = cond_B_design, SL1344 = cond_C_design, IFNg_SL1344 = cond_D_design)
  
#### EXPRESSION ####
#Filter expression data by samples
exprs_cqn = expression_dataset$exprs_cqn[,design$sample_id]
exprs_counts = expression_dataset$exprs_counts[,design$sample_id]

#Define expressed genes
mean_expression = calculateMean(exprs_cqn, as.data.frame(design), "condition_name")
expressed_genes = names(which(apply(mean_expression, 1, max) > 0))

#Filter expression data by min expression
exprs_cqn = expression_dataset$exprs_cqn[expressed_genes, design$sample_id]
exprs_counts = expression_dataset$exprs_counts[expressed_genes, design$sample_id]

#Set up expression data for each condition
condA_exp = extractSubset(dplyr::filter(sample_meta, condition == "A"), exprs_cqn)
condB_exp = extractSubset(dplyr::filter(sample_meta, condition == "B"), exprs_cqn)
condC_exp = extractSubset(dplyr::filter(sample_meta, condition == "C"), exprs_cqn)
condD_exp = extractSubset(dplyr::filter(sample_meta, condition == "D"), exprs_cqn)
exprs_cqn_list = list(naive = condA_exp, IFNg = condB_exp, SL1344 = condC_exp, IFNg_SL1344 = condD_exp)

#Set up the genepos data.frame
genepos = dplyr::filter(expression_dataset$gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::transmute(chr = chromosome_name, left = start_position, right = end_position, geneid = gene_id, score = 1000, strand) %>%
  dplyr::mutate(strand = ifelse(strand == 1, "+","-")) %>%
  as.data.frame()

#### GENOTYPES ####
#Remove two SNPs that have the same name for two different positions
duplicates_count = table(vcf_file$snpspos$snpid)
duplicate_ids = names(which(duplicates_count > 1))
vcf_file$snpspos = dplyr::filter(vcf_file$snpspos, !(snpid %in% duplicate_ids))
vcf_file$genotypes = vcf_file$genotypes[!(rownames(vcf_file$genotypes) %in% duplicate_ids),]

#Filter genotype data
geno_data = extractSubset(dplyr::filter(sample_meta, condition == "A"), vcf_file$genotypes, old_column_names = "genotype_id")

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
  genotypes = geno_data,
  snpspos = vcf_file$snpspos,
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

