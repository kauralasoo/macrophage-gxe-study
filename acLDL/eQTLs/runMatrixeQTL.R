library("MatrixEQTL")
library("ggplot2")
library("devtools")
load_all("macrophage-gxe-study/seqUtils/")

#Load expression dataset
dataset = readRDS("results/acLDL/acLDL_combined_expression_data.rds")
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/acLDL/acLDL_array_genotypes.GRCh38.vcfToMatrix.rds") #genotypes

#Remove outlier samples xegx and nusw
filtered_design = dplyr::filter(dataset$design, !(donor %in% c("nusw","xegx")))
rownames(filtered_design) = filtered_design$sample_id
filtered_exprs_cqn = dataset$exprs_cqn[,filtered_design$sample_id]
filtered_counts = dataset$exprs_counts[,filtered_design$sample_id]

#Extraxt donor genotype match from metadata
donor_genotype_match = dplyr::left_join(filtered_design, line_metadata, by = "donor") %>% 
  dplyr::select(line_id, donor, genotype_id) %>% unique()
sample_meta = dplyr::left_join(filtered_designc, donor_genotype_match, by = "donor")

#Filter genes by expression
expressed_genes = names(which(rowMeans(filtered_exprs_cqn) > 0)) #Set conservative threshold to expression level
exprs_sel = filtered_exprs_cqn[expressed_genes,]

#Correct vorx and zuta swap
indexes = which(colnames(exprs_sel) == "ZUTA_24h_Ctrl" | colnames(exprs_sel) == "ZUTA_24h_AcLDL" | colnames(exprs_sel) == "VORX_24h_Ctrl" | colnames(exprs_sel) == "VORX_24h_AcLDL")
colnames(exprs_sel)[indexes] = c("ZUTA_24h_AcLDL","ZUTA_24h_Ctrl","VORX_24h_AcLDL","VORX_24h_Ctrl")

#Filter genotype data
sample_meta_filtered = dplyr::filter(sample_meta, condition == "Ctrl")
geno_data = qtlProcessGenotypes(sample_meta_filtered, vcf_file$genotypes)

#Filter expression data
ctrl_exp = qtlProcessExpression(dplyr::filter(sample_meta, condition == "Ctrl"), exprs_sel)
acldl_exp = qtlProcessExpression(dplyr::filter(sample_meta, condition == "AcLDL"), exprs_sel)

#Set up the genepos data.frame
genepos = dplyr::filter(dataset$gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::transmute(geneid = gene_id, chr = chromosome_name, left = start_position, right = end_position) %>%
  as.data.frame()

#Run MatrixeQTL
ctrl_res = runMatrixEQTL(ctrl_exp, geno_data, vcf_file$snpspos, genepos)
write.table(ctrl_res$cis$eqtls, "results/acLDL/matrixeQTL/matrixeQTL_Ctrl.txt", sep ="\t", quote = FALSE, row.names = FALSE)
acldl_res = runMatrixEQTL(acldl_exp, geno_data, vcf_file$snpspos, genepos)
write.table(acldl_res$cis$eqtls, "results/acLDL/matrixeQTL/matrixeQTL_AcLDL.txt", sep ="\t", quote = FALSE, row.names = FALSE)

