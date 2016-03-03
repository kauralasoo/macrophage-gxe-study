library("devtools")
load_all("../seqUtils/")
library("dplyr")
library("ggplot2")
library("tidyr")

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")
snp_positions = vcf_file$snpspos %>% tbl_df() %>% dplyr::rename(snp_id = snpid) 
saveRDS(snp_positions, "genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.snp_positions.rds")
snp_positions = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.snp_positions.rds")

gwas_catalog = importGwasCatalog("annotations/gwas_catalog_v1.0.1-downloaded_2016-03-02.tsv") %>% 
  dplyr::rename(catalog_snp_id = snp_id)

#Find matching SNPs
match_by_coords = dplyr::semi_join(snp_positions, gwas_catalog, by = c("chr","pos"))

#Filter by matching SNPs
matched_catalog = dplyr::left_join(gwas_catalog, match_by_coords, by = c("chr","pos")) %>%
  dplyr::filter(!is.na(snp_id)) %>% #Has a matching SNP in our genotype data
  dplyr::filter(is_european == TRUE) %>% #Only keep European studies
  dplyr::filter(sample_size > 1000) %>% #At least 1000 samples
  dplyr::filter(gwas_pvalue < 1) %>% #Obvious p-value errors
  dplyr::group_by(trait, snp_id) %>% 
  dplyr::arrange(gwas_pvalue) %>%
  dplyr::filter(row_number() == 1) %>% #Keep single copy of a SNP per trair
  dplyr::ungroup()

#Calculate R2 between the SNPs
trait_snp_pairs = dplyr::transmute(matched_catalog, gene_id = trait, snp_id)

#Filter the VCF file
filtered_genotypes = vcf_file$genotypes[unique(trait_snp_pairs$snp_id),]
filtered_pairs = filterHitsR2(trait_snp_pairs, filtered_genotypes, .8) %>%
  dplyr::rename(trait = gene_id)
filtered_catalog = dplyr::semi_join(matched_catalog, filtered_pairs, by = c("trait","snp_id"))

#Save catalog to disk
saveRDS(filtered_catalog, "macrophage-gxe-study/data/gwas_catalog/gwas_catalog_v1.0.1-downloaded_2016-03-02.filtered.rds")



