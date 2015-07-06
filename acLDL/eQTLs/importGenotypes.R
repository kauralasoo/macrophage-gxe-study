library("dplyr")
library("tidyr")
library("devtools")
library("VariantAnnotation")
load_all("macrophage-gxe-study/seqUtils/")

#Import genotype data from the VCF file
vcf_file = vcfToMatrix("genotypes/acLDL/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted.filtered.vcf.gz", "GRCh38")
saveRDS(vcf_file, "genotypes/acLDL/acLDL_array_genotypes.GRCh38.vcfToMatrix.rds")

