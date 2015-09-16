library("MatrixEQTL")
library("ggplot2")
library("devtools")
library("dplyr")
load_all("../seqUtils/")

#Load expresison dataset preapred previously by processExpressionData.R script
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/SL1344/array_genotypes.59_samples.vcfToMatrix.rds") #genotypes

#Discard replicate samples
design = dplyr::filter(expression_dataset$design, !(donor == "fpdj")) %>% tbl_df() %>% #Remove all fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) %>% #Remove second ougl sample (dium)
  dplyr::filter(!(donor == "mijn")) %>% #Remove mijn (wrong line from CGAP)
  dplyr::filter(!(donor == "jorr")) #Very strong outlier in PEER analysis

sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))

#Export data for PEER
exprs_cqn = expression_dataset$exprs_cqn[,design$sample_id]
expressed_genes = names(which(rowMeans(exprs_cqn) > 0)) #Set conservative threshold to expression level
exprs_cqn = expression_dataset$exprs_cqn[expressed_genes,design$sample_id]

#Condition A
cond_A_design = dplyr::filter(design, condition == "A")
cond_A_exprs = t(exprs_cqn[,cond_A_design$sample_id])
write.table(cond_A_exprs, "results/SL1344/cond_A_exprs.peer.txt", row.names = FALSE, col.names = FALSE, sep = ",")

#Filter expression data by min expression
expressed_genes = names(which(rowMeans(exprs_cqn) > 0)) #Set conservative threshold to expression level
exprs_sel = exprs_cqn[expressed_genes,]

#Filter genotype data
sample_meta_filtered = dplyr::filter(sample_meta, condition == "A")
geno_data = qtlProcessGenotypes(sample_meta_filtered, vcf_file$genotypes)

#Set up expression data for each condition
condA_exp = qtlProcessExpression(dplyr::filter(sample_meta, condition == "A"), exprs_sel)
condB_exp = qtlProcessExpression(dplyr::filter(sample_meta, condition == "B"), exprs_sel)
condC_exp = qtlProcessExpression(dplyr::filter(sample_meta, condition == "C"), exprs_sel)
condD_exp = qtlProcessExpression(dplyr::filter(sample_meta, condition == "D"), exprs_sel)

#Set up the genepos data.frame
genepos = dplyr::filter(expression_dataset$gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::transmute(geneid = gene_id, chr = chromosome_name, left = start_position, right = end_position) %>%
  as.data.frame()

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

