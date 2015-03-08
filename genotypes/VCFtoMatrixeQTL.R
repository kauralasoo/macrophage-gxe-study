#Convert a VCF file into a format suitable for MatrixEQTL package
library("VariantAnnotation")
library("dplyr")

#Load the VCF file from disk
genotypes_vcf = readVcf("genotypes/selected_genotypes.GRCh38.sorted.vcf.gz", "GRCh38")

# Extract SNP positions from the VCF file
variant_granges = rowData(genotypes_vcf)
elementMetadata(variant_granges) = c()
snp_positions = as.data.frame(variant_granges)
snp_df = dplyr::mutate(snp_positions, snpid = rownames(snp_positions)) %>% 
  tbl_df() %>%
  dplyr::select(snpid, seqnames, start) %>%
  dplyr::rename(chr = seqnames, pos = start)
write.table(snp_df, "genotypes/snp_positions.tsv", sep ="\t", quote = FALSE, row.names = FALSE)
saveRDS(snp_df, "genotypes/snp_positions.rds")

#Extract genotype matrix
genotypes = geno(genotypes_vcf)$GT
genotypes[genotypes == "1/1"] = 2
genotypes[genotypes == "0/1"] = 1
genotypes[genotypes == "1/0"] = 1
genotypes[genotypes == "0/0"] = 0
genotypes[genotypes == "."] = "NA"
write.table(genotypes, "genotypes/genotypes_for_matrixeQTL.tsv", sep ="\t", quote = FALSE)
saveRDS(genotypes, "genotypes/genotypes_for_matrixeQTL.rds")
