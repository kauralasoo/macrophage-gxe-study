load_all("macrophage-gxe-study/seqUtils/")
library(plyr)
library(dplyr)

#Load the raw eQTL dataset
eqtl_dataset = readRDS("results/SL1344/matrixeQTL/matrixeQTL_dataset.rds")

#Load identified eQTLS
condA_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_A.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE))
condB_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_B.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE))
condC_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_C.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE))
condD_eqtls = filterEQTLs(read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_D.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE))


#eQTLS found in SL1344
condA_snps = dplyr::select(condA_eqtls, gene_id = gene, genotype_id = snps)
condA_plots = makeMultiplePlots(condA_snps, eqtl_dataset)
savePlots(condA_plots, path = "results/SL1344/matrixeQTL/condA_plots/",width = 7, height = 7)

#eQTLS found in IFNg
condB_snps = dplyr::select(condB_eqtls, gene_id = gene, genotype_id = snps)
condB_plots = makeMultiplePlots(condB_snps, eqtl_dataset)
savePlots(condB_plots, path = "results/SL1344/matrixeQTL/condB_plots/",width = 7, height = 7)

#eQTLS found in SL1344
condC_snps = dplyr::select(condC_eqtls, gene_id = gene, genotype_id = snps)
condC_plots = makeMultiplePlots(condC_snps, eqtl_dataset)
savePlots(condC_plots, path = "results/SL1344/matrixeQTL/condC_plots/",width = 7, height = 7)

#eQTLS found in SL1344
condD_snps = dplyr::select(condD_eqtls, gene_id = gene, genotype_id = snps)
condD_plots = makeMultiplePlots(condD_snps, eqtl_dataset)
savePlots(condD_plots, path = "results/SL1344/matrixeQTL/condD_plots/",width = 7, height = 7)


condA_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_A.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
condB_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_B.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
condC_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_C.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
condD_eqtls = read.table("results/SL1344/matrixeQTL/temp/matrixeQTL_condition_D.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
eqtl_list = list(A = condA_eqtls, B = condB_eqtls, C = condC_eqtls, D = condD_eqtls)
