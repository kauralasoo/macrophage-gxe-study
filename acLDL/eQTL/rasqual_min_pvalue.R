library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("Rsamtools")

#Import data
acldl_list = readRDS("results/acLDL/acLDL_combined_expression_data_covariates.rds")


#Import SNP coordinates
snp_coords = readr::read_delim("genotypes/acLDL/imputed_20151005/imputed.70_samples.snp_coords.txt", 
                               delim = "\t", col_types = "cdc", col_names = c("chr","pos","snp_id"))

#Extract minimal p-value for each condition
ctrl_eigen_pvalue = eigenMTImportResults("results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.eigenMT.txt")
acldl_eigen_pvalue = eigenMTImportResults("results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.eigenMT.txt")

min_pvalue_list = list(Ctrl = ctrl_eigen_pvalue,
                       AcLDL = acldl_eigen_pvalue)
saveRDS(min_pvalue_list, "results/acLDL/eQTLs/acLDL_rasqual_min_pvalues.rds")