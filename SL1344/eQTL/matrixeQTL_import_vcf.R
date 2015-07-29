library("dplyr")
library("tidyr")
library("devtools")
library("cqn")
load_all("macrophage-gxe-study/seqUtils/")

#Import genotype data from the VCF file
vcf_file = vcfToMatrix("genotypes/SL1344/array_genotypes.59_samples.vcf.gz", "GRCh38")
saveRDS(vcf_file, "genotypes/SL1344/array_genotypes.59_samples.vcfToMatrix.rds")

