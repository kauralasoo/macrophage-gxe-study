library("devtools")
library("SNPRelate")
library("GWASTools")
load_all("../seqUtils/")

#Import VCF file into R and convert it into a matrix of allele counts
#All variants
#SNPRelate::snpgdsVCF2GDS("processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz", 
#                         "processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.gds", method = "copy.num.of.ref")
vcf_file = seqUtils::gdsToMatrix("processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.gds")
saveRDS(vcf_file, "processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.rds")
