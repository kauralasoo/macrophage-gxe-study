library("dplyr")

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#IBD summary stats
ibd_hits = readr::read_delim("databases/GWAS/IBD/IBD_trans_ethnic_association_summ_stats_b37.txt", delim = " ")

#Filter IBD hits
ibd_filtered = dplyr::transmute(ibd_hits, ibd_snp_id = SNP, A1_effect, A2_other,beta_EUR, se_EUR, P_EUR)

#import new coords
ibd_new_coords = readr::read_delim("databases/GWAS/IBD/IBD_trans_ethnic_coords.GRCh38.bed", delim = "\t", 
                                   col_names = c("chr","pos", "pos2", "ibd_snp_id")) %>%
  dplyr::select(chr, pos, ibd_snp_id)

ibd_matched = dplyr::left_join(ibd_new_coords, ibd_filtered, by = "ibd_snp_id") %>%
  dplyr::mutate(chr = as.character(chr))

#Extract SNP ids from VCF file
new_snp_ids = dplyr::semi_join(vcf_file$snpspos, ibd_matched, by = c("chr", "pos")) %>%
  dplyr::rename(snp_id = snpid)

#Add current SNP ids
ibd_final = dplyr::left_join(ibd_matched, new_snp_ids, by = c("chr","pos"))
saveRDS(ibd_final, "annotations/IBD_trans_ethnic_summary_stats.rds")

#Export stats in GWAS format
ibd_gwas = dplyr::transmute(ibd_final, CHR = chr, BP = pos, SNP = ibd_snp_id, P = P_EUR) %>%
  dplyr::arrange(CHR, BP)
write.table(ibd_gwas, "annotations/IBD_trans_ethnic_summary_stats.gwas", sep ="\t", quote = FALSE, row.names = FALSE)

