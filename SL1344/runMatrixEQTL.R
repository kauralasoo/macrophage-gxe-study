library("MatrixEQTL")
library("ggplot2")
load_all("macrophage-gxe-study/seqUtils/")

#Load expresison dataset preapred previously by processExpressionData.R script
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/selected_genotypes.GRCh38.vcfToMatrix.rds") #genotypes

#Discard replicate samples
design = dplyr::filter(expression_dataset$design, !(donor == "ffdp")) %>% #ffdp identical to fpdl
  dplyr::filter(!(donor == "fpdj" & replicate == 2)) #fpdj has two replicates
sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate"))

#Filter expression data by min expression
exprs_cqn = expression_dataset$exprs_cqn[,design$sample_id]
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
condD_res = runMatrixEQTL(condD_exp, geno_data, vcf_file$snpspos, genepos)write.table(condD_res$cis$eqtls, "results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_D.txt", sep ="\t", quote = FALSE, row.names = FALSE)
write.table(condD_res$cis$eqtls, "results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_D.txt", sep ="\t", quote = FALSE, row.names = FALSE)


