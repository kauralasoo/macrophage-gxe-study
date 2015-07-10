load_all("macrophage-gxe-study/seqUtils/")
library(plyr)
library(dplyr)

#Load the raw eQTL dataset
expression_dataset = readRDS("results/SL1344/combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/selected_genotypes.GRCh38.vcfToMatrix.rds") #genotypes
gene_id_name_map = dplyr::select(expression_dataset$gene_metadata, gene_id, gene_name)

#Load identified eQTLS
condA_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_A.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map)
condB_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_B.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map)
condC_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_C.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map)
condD_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/fixed/matrixeQTL_condition_D.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map)


#eQTLS found in SL1344
condA_plots = makeMultiplePlots(condA_eqtls, expression_dataset, vcf_file, line_metadata)
savePlots(condA_plots, path = "results/SL1344/matrixeQTL/condA_plots/",width = 7, height = 7)

#eQTLS found in IFNg
condB_plots = makeMultiplePlots(condB_eqtls, expression_dataset, vcf_file, line_metadata)
savePlots(condB_plots, path = "results/SL1344/matrixeQTL/condB_plots/",width = 7, height = 7)

#eQTLS found in SL1344
condC_plots = makeMultiplePlots(condC_eqtls, expression_dataset, vcf_file, line_metadata)
savePlots(condC_plots, path = "results/SL1344/matrixeQTL/condC_plots/",width = 7, height = 7)

#eQTLS found in SL1344
condD_plots = makeMultiplePlots(condD_eqtls, expression_dataset, vcf_file, line_metadata)
savePlots(condD_plots, path = "results/SL1344/matrixeQTL/condD_plots/",width = 7, height = 7)


condA_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_A.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
condB_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_B.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
condC_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_C.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
condD_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_D.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
eqtl_list = list(A = condA_eqtls, B = condB_eqtls, C = condC_eqtls, D = condD_eqtls)
