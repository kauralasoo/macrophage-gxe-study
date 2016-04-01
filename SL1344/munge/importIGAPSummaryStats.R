#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#STAGE 1_2
#Import IGAP stats
igap_raw = read.table("databases/GWAS/IGAP_summary_statistics/IGAP_stage_1_2_combined.txt", header = TRUE, stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  dplyr::rename(igap_snp_id = MarkerName, igap_pvalue = Pvalue)
#Import new coords
igap_new_coords = read.table("databases/GWAS/IGAP_summary_statistics/IGAP_stage_1_2_coords.GRCh38.bed", stringsAsFactors = FALSE)  %>% 
  tbl_df() %>% dplyr::transmute(chr = V1, pos = V2, igap_snp_id = V4)

#Line the two together
igap_processed = dplyr::left_join(igap_new_coords, igap_raw, by = "igap_snp_id") %>% 
  dplyr::select(chr, pos, igap_snp_id, Effect_allele, Non_Effect_allele, Beta, SE, igap_pvalue) %>%
  dplyr::mutate(chr = as.character(chr))

#Extract SNP ids from VCF file
new_snp_ids = dplyr::semi_join(vcf_file$snpspos, igap_processed, by = c("chr", "pos")) %>%
  dplyr::rename(snp_id = snpid)

#Add current SNP ids
igap_final = dplyr::left_join(igap_processed, new_snp_ids, by = c("chr","pos"))
saveRDS(igap_final, "annotations/IGAP_stage_1_2_combined_stats.rds")

#Export stats in GWAS format
igap_gwas = dplyr::transmute(igap_final, CHR = chr, BP = pos, SNP = igap_snp_id, P = igap_pvalue)
write.table(igap_gwas, "annotations/IGAP_stage_1_2_combined_stats.gwas", sep ="\t", quote = FALSE, row.names = FALSE)

#STAGE 1 only
#Import IGAP stats
igap_raw = readr::read_delim("databases/GWAS/IGAP_summary_statistics/IGAP_stage_1.txt", delim = "\t") %>%
  tbl_df() %>%
  dplyr::rename(igap_snp_id = MarkerName, igap_pvalue = Pvalue)
#Import new coords
igap_new_coords = readr::read_delim("databases/GWAS/IGAP_summary_statistics/IGAP_stage_1_coords.GRCh38.bed", delim = "\t", col_names = FALSE) %>% 
  dplyr::transmute(chr = X1, pos = X2, igap_snp_id = X4)

#Line the two together
igap_processed = dplyr::left_join(igap_new_coords, igap_raw, by = "igap_snp_id") %>% 
  dplyr::select(chr, pos, igap_snp_id, Effect_allele, Non_Effect_allele, Beta, SE, igap_pvalue) %>%
  dplyr::mutate(chr = as.character(chr))

#Line the two together
igap_processed = dplyr::left_join(igap_new_coords, igap_raw, by = "igap_snp_id") %>% 
  dplyr::select(chr, pos, igap_snp_id, Effect_allele, Non_Effect_allele, Beta, SE, igap_pvalue) %>%
  dplyr::mutate(chr = as.character(chr))

