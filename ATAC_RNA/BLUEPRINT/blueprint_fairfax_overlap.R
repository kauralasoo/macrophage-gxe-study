library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")

importBlueprintSummary <- function(summary_path){
  colnames = c("snp_string", "old_snp_id", "peak_id","p_nominal","beta", "p_bonferroni", "FDR", "alt_AF","std_error")
  data = readr::read_delim(summary_path, delim = "\t", col_names = colnames, col_types = "cccdddddd") %>%
    dplyr::filter(p_bonferroni <= 1) %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, "fdr")) %>%
    tidyr::separate(snp_string, c("chr","pos_string"), sep = ":") %>%
    tidyr::separate(pos_string, c("pos","ref","alt"), sep = "_") %>%
    dplyr::mutate(pos = as.integer(pos))
  return(data)
}

#Import BLUEPRINT min p-values
k4me1 = importBlueprintSummary("databases/BLUEPRINT/mono_K4ME1_min_pvalues.txt.gz")
k27ac = importBlueprintSummary("databases/BLUEPRINT/mono_K27AC_min_pvalues.txt.gz")

#Import GRCh37 variant information
variant_information = importVariantInformation("processed/Fairfax/merged_genotypes/fairfax_genotypes.variant_information.txt.gz")

#Filter k4me1 and k27ac results for SNPs that have genotype data
new_snps = dplyr::select(variant_information, chr, pos, snp_id) %>% dplyr::filter(pos %in% k4me1$pos)
k4me1_filtered = dplyr::left_join(k4me1, new_snps, by = c("chr","pos")) %>% dplyr::filter(!is.na(snp_id))

new_snps = dplyr::select(variant_information, chr, pos, snp_id) %>% dplyr::filter(pos %in% k27ac$pos)
k27ac_filtered = dplyr::left_join(k27ac, new_snps, by = c("chr","pos")) %>% dplyr::filter(!is.na(snp_id))

#Construct metadata for chromatin qtls
k4me1_meta = dplyr::select(k4me1_filtered, peak_id) %>% 
  tidyr::separate(peak_id, c("chr","start","end"), sep = ":", remove = FALSE, convert = TRUE)
k27ac_meta = dplyr::select(k27ac_filtered, peak_id) %>% 
  tidyr::separate(peak_id, c("chr","start","end"), sep = ":", remove = FALSE, convert = TRUE)

#Filter QTLs
k4me1_qtls = dplyr::filter(k4me1_filtered, p_fdr < 0.1)
k27ac_qtls = dplyr::filter(k27ac_filtered, p_fdr < 0.1)

#Import Fairfax QTLs
#Import SummarizedExperiment
se_fairfax = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")
gene_name_map = rowData(se_fairfax) %>% tbl_df2() %>% dplyr::mutate(gene_id = probe_id) %>%
  dplyr::select(gene_id, gene_name, chr)

#Import aFC estimates
shared_diff_df = readRDS("results/Fairfax/interactions/shared_aFC_estimates.rds") %>%
  dplyr::filter(abs(CD14) > 0.29)

#Import QTLs
all_qtls = readRDS("results/Fairfax/fairfax_qtl_min_pvalues.rds")$shared$CD14 %>%
  dplyr::filter(p_fdr < 0.1) %>%
  dplyr::rename(gene_id = phenotype_id) %>%
  dplyr::semi_join(shared_diff_df, by = c("gene_id", "snp_id"))
response_qtls = readRDS("results/Fairfax/interactions/shared_interaction_hits.rds") %>%
  dplyr::left_join(gene_name_map, by = c("gene_id", "gene_name"))

#Find overlaps chr-by-chr
chr_list = c(1:22) %>% as.character() %>% idVectorToList() 

#Find overlaps with response qtls
response_olap_list = purrr::map(chr_list, function(chrom){
  print(chrom)
  vcf_file = seqUtils::gdsToMatrix(paste0("processed/Fairfax/geno_by_chr/",chrom,".gds"))
  selected_qtls = dplyr::filter(response_qtls, chr == chrom)
  
  k4me1_overlaps = findGWASOverlaps(dplyr::select(selected_qtls, gene_id, snp_id), 
                                    dplyr::select(k4me1_qtls, peak_id, snp_id, chr, pos), 
                                    vcf_file, max_distance = 5e5, min_r2 = 0.8)
  
  k27ac_overlaps = findGWASOverlaps(dplyr::select(selected_qtls, gene_id, snp_id), 
                                    dplyr::select(k27ac_qtls, peak_id, snp_id, chr, pos), 
                                    vcf_file, max_distance = 5e5, min_r2 = 0.8)
  return(list(k4me = k4me1_overlaps, k27ac = k27ac_overlaps))
})

k27ac_overlaps = purrr::map_df(response_olap_list, ~.$k27ac)
k4me_overlaps = purrr::map_df(response_olap_list, ~.$k4me)

length(unique(k27ac_overlaps$gene_id)) / length(unique(response_qtls$gene_id))
length(unique(k4me_overlaps$gene_id)) / length(unique(response_qtls$gene_id))

length(unique(k27ac_overlaps$gene_id))/length(unique(c(k27ac_overlaps$gene_id,k4me_overlaps$gene_id)))


#Find overlaps with all eQTLs
all_olap_list = purrr::map(chr_list, function(chrom){
  print(chrom)
  vcf_file = seqUtils::gdsToMatrix(paste0("processed/Fairfax/geno_by_chr/",chrom,".gds"))
  selected_qtls = dplyr::filter(all_qtls, pheno_chr == chrom)
  
  k4me1_overlaps = findGWASOverlaps(dplyr::select(selected_qtls, gene_id, snp_id), 
                                    dplyr::select(k4me1_qtls, peak_id, snp_id, chr, pos), 
                                    vcf_file, max_distance = 5e5, min_r2 = 0.8)
  
  k27ac_overlaps = findGWASOverlaps(dplyr::select(selected_qtls, gene_id, snp_id), 
                                    dplyr::select(k27ac_qtls, peak_id, snp_id, chr, pos), 
                                    vcf_file, max_distance = 5e5, min_r2 = 0.8)
  return(list(k4me = k4me1_overlaps, k27ac = k27ac_overlaps))
})

k27ac_olaps_all = purrr::map_df(all_olap_list, ~.$k27ac)
k4me_olaps_all = purrr::map_df(all_olap_list, ~.$k4me)

length(unique(k27ac_olaps_all$gene_id)) / length(unique(all_qtls$gene_id))
length(unique(k4me_olaps_all$gene_id)) / length(unique(all_qtls$gene_id))

length(unique(k27ac_olaps_all$gene_id))/length(unique(c(k27ac_olaps_all$gene_id,k4me_olaps_all$gene_id)))

fisher.test(matrix(c(152, 226-152,378, 469-378), ncol = 2))


