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

#Set up expression data for each condition
condA_exp = extractSubset(dplyr::filter(sample_meta, condition == "A"), exprs_cqn)
condB_exp = extractSubset(dplyr::filter(sample_meta, condition == "B"), exprs_cqn)
condC_exp = extractSubset(dplyr::filter(sample_meta, condition == "C"), exprs_cqn)
condD_exp = extractSubset(dplyr::filter(sample_meta, condition == "D"), exprs_cqn)

#Make list of expression levels
expression_list = list(naive = condA_exp, IFNg = condB_exp, SL1344 = condC_exp, IFNg_SL1344 = condD_exp)
saveEqtlbmaData(expression_list, output_dir = "eqtlbma/input", "expression")

#Construct covariate matrix
covariates = dplyr::mutate(sample_meta, sex = ifelse(gender == "male",1,0)) %>%
  dplyr::select(sample_id, sex, ng_ul_mean, diff_days)

#Import PEER factors
peer_factors = readRDS("results/SL1344/PEER/expressed_10_factors.rds")
peer_covariates = dplyr::left_join(covariates, peer_factors, by = "sample_id") %>% as.data.frame()
rownames(peer_covariates) = peer_covariates$sample_id
peer_covariates = t(peer_covariates[,-1])

#Set up covatiates for each conditon
condA_covs = extractSubset(dplyr::filter(sample_meta, condition == "A"), peer_covariates)[1:7,]
condB_covs = extractSubset(dplyr::filter(sample_meta, condition == "B"), peer_covariates)[1:7,]
condC_covs = extractSubset(dplyr::filter(sample_meta, condition == "C"), peer_covariates)[1:7,]
condD_covs = extractSubset(dplyr::filter(sample_meta, condition == "D"), peer_covariates)[1:7,]

#Make list of covariates
covariates_list = list(naive = condA_covs, IFNg = condB_covs, SL1344 = condC_covs, IFNg_SL1344 = condD_covs)
saveEqtlbmaData(covariates_list, output_dir = "eqtlbma/input", "covariates")

#Genotypes
#Filter genotype data
sample_meta_filtered = dplyr::filter(sample_meta, condition == "A")
geno_data = extractSubset(sample_meta_filtered, vcf_file$genotypes, old_column_names = "genotype_id")
write.table(geno_data, "eqtlbma/input/genotypes.txt", sep ="\t", quote = FALSE)

#Set up the genepos data.frame
genepos = dplyr::filter(expression_dataset$gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::transmute(geneid = gene_id, chr = chromosome_name, left = start_position, right = end_position) %>%
  as.data.frame()

#Gene positions
genepos_eqtlbma = dplyr::filter(expression_dataset$gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::transmute(chr = chromosome_name, left = start_position, right = end_position, geneid = gene_id, score = 1000, strand) %>%
  dplyr::mutate(strand = ifelse(strand == 1, "+","-")) %>%
  as.data.frame()
write.table(genepos_eqtlbma, "eqtlbma/input/gene_coords.bed", sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#SNP positions
snpspos_eqtlbma = dplyr::transmute(vcf_file$snpspos, chr, start  = pos-1, end = pos, snpid)
write.table(snpspos_eqtlbma, "eqtlbma/input/snp_coords.bed", sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#Filter genotype data
sample_meta_filtered = dplyr::filter(sample_meta, condition == "A")
geno_data = extractSubset(sample_meta_filtered, vcf_file$genotypes, old_column_names = "genotype_id")

#Set up the genepos data.frame
genepos = dplyr::filter(expression_dataset$gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::transmute(geneid = gene_id, chr = chromosome_name, left = start_position, right = end_position) %>%
  as.data.frame()


