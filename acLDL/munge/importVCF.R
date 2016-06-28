library("devtools")
library("SNPRelate")
load_all("../seqUtils/")

#Import VCF file into R and convert it into a matrix of allele counts
#All variants
SNPRelate::snpgdsVCF2GDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.vcf.gz", 
                         "genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.gds", method = "copy.num.of.ref")
vcf_file = seqUtils::gdsToMatrix("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.gds")
saveRDS(vcf_file, "genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.rds")

#Variants with INFO > 0.7
SNPRelate::snpgdsVCF2GDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz", 
                         "genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.gds", method = "copy.num.of.ref")
vcf_file = seqUtils::gdsToMatrix("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.gds")
saveRDS(vcf_file, "genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.rds")
