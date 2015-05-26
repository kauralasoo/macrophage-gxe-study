library("dplyr")
library("tidyr")
library("devtools")
library("cqn")
load_all("macrophage-gxe-study/seqUtils/")

#Import genotype data from the VCF file
vcf_file = vcfToMatrix("genotypes/selected_genotypes.GRCh38.sorted.vcf.gz", "GRCh38")
saveRDS(vcf_file, "genotypes/selected_genotypes.GRCh38.vcfToMatrix.rds")

