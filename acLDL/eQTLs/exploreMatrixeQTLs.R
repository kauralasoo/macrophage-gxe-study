library("ggplot2")
library("devtools")
load_all("macrophage-gxe-study/seqUtils/")

#Load expression dataset
dataset = readRDS("results/acLDL/acLDL_combined_expression_data.rds")
colnames(dataset$design)[4] = "condition_name"
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
vcf_file = readRDS("genotypes/acLDL/acLDL_array_genotypes.GRCh38.vcfToMatrix.rds") #genotypes
gene_id_name_map = dplyr::select(dataset$gene_metadata, gene_id, gene_name)

#Load identified eQTLS
ctrl_eqtls = filterEQTLs(read.table("results/acLDL/matrixeQTL/matrixeQTL_Ctrl.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map)
acldl_eqtls = filterEQTLs(read.table("results/acLDL/matrixeQTL/matrixeQTL_AcLDL.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE), gene_id_name_map)

#eQTLS found in SL1344
ctrl_plots = makeMultiplePlots_acLDL(ctrl_eqtls, dataset, vcf_file, line_metadata)
savePlots(ctrl_plots, path = "results/acLDL/matrixeQTL/ctrl_plots/",width = 7, height = 7)

#eQTLS found in SL1344
acldl_plots = makeMultiplePlots_acLDL(acldl_eqtls, dataset, vcf_file, line_metadata)
savePlots(acldl_plots, path = "results/acLDL/matrixeQTL/acldl_plots/",width = 7, height = 7)
