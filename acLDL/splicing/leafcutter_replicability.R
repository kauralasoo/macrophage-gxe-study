library("plyr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("Rsamtools")
library("purrr")

#Import data
prop_list = readRDS("results/acLDL/acLDL_combined_proportions.row_quantile.rds")

#Import variant information
snp_info = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

#Find minimal p-values from fastQTL results
ctrl_fqtl = importFastQTLTable("results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_permuted.txt.gz")
acldl_fqtl = importFastQTLTable("results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_permuted.txt.gz")

fastqtl_pvalue_list = list(Ctrl = ctrl_fqtl,
                           AcLDL = acldl_fqtl)
saveRDS(fastqtl_pvalue_list, "results/acLDL/leafcutter/leafcutter_min_pvalues.rds")
