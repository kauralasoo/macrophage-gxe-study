library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import GRCh37 variant information
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")

#Import eQTLs
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")

#Iterate over chromosomes
chromosomes = c(1:22) %>% as.character()
for (chr in chromosomes){
  chr_string = paste0("chr", chr)
  print(chr_string)
  
  #Import Garfield coords
  variant_coords_dir = "databases/garfield-data/annotations/variant_coords/"
  coords = readr::read_tsv(file.path(variant_coords_dir, chr_string), col_names = "pos", col_types = "i")
  
  #Process summary stats
  garfield_results = rasqualSummariesToGarfieldByChr(chr, rasqual_min_pvalues, qtlResults()$rna_rasqual, 
                        GRCh37_variants, combined_expression_data$gene_metadata, coords)
  
  #Export annotations
  output_dir = "databases/garfield-data/annotations/"
  write.table(garfield_results, file.path(output_dir, chr_string), sep = "", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  
}

#Import caQTLs
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
rasqual_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")

#Iterate over chromosomes
chromosomes = c(1:22) %>% as.character()
for (chr in chromosomes){
  chr_string = paste0("chr", chr)
  print(chr_string)
  
  #Import Garfield coords
  variant_coords_dir = "databases/garfield-data/annotations/variant_coords/"
  coords = readr::read_tsv(file.path(variant_coords_dir, chr_string), col_names = "pos", col_types = "i")
  
  #Process summary stats
  garfield_results = rasqualSummariesToGarfieldByChr(chr, rasqual_min_pvalues, qtlResults()$atac_rasqual, 
                                                     GRCh37_variants, atac_data$gene_metadata, coords)
  
  #Export annotations
  output_dir = "databases/garfield-data/annotation_caqtl/"
  write.table(garfield_results, file.path(output_dir, chr_string), sep = "", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  
}




#Export only lead eQTL and caQTL variants for GARFIELD
eQTL_rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
names(eQTL_rasqual_min_pvalues) = paste0(names(eQTL_rasqual_min_pvalues), "_eQTL")

caQTL_rasqual_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
names(caQTL_rasqual_min_pvalues) = paste0(names(caQTL_rasqual_min_pvalues), "_caQTL")

#Merge all QTLs
all_qtls_snps = c(eQTL_rasqual_min_pvalues, caQTL_rasqual_min_pvalues) %>% 
  purrr::map(~dplyr::filter(., p_fdr < 0.1) %>% dplyr::select(snp_id) %>% unique())

#Add GRCh37 coordinates
lead_variants = purrr::map_df(all_qtls_snps, identity) %>% unique()
selected_variants = dplyr::filter(GRCh37_variants, snp_id %in% lead_variants$snp_id) %>% 
  dplyr::select(chr, pos, snp_id)

#Add variants
all_qtl_pos = purrr::map(all_qtls_snps, ~dplyr::left_join(., selected_variants, by = "snp_id") %>%
                           dplyr::filter(!is.na(pos)) %>%
                           dplyr::mutate(is_qtl = 1))

#Extract all variants
all_variants = purrr::map_df(all_qtl_pos, identity) %>% unique() %>%
  dplyr::select(snp_id, chr, pos)

#Merge all lists into a single df
reduced_df = purrr::reduce(all_qtl_pos, dplyr::left_join, .init = all_variants, by = c("snp_id", "chr", "pos"))
names(reduced_df)[4:11] = names(all_qtl_pos)
reduced_df[is.na(reduced_df)] = 0

#Iterate over chromosomes and save results
#Iterate over chromosomes
chromosomes = c(1:22) %>% as.character()
for (chr in chromosomes){
  chr_string = paste0("chr", chr)
  print(chr_string)
  
  #Extract chromosome data 
  chr_df = extractChrFromReducedDf(chr, reduced_df)
  
  output_dir = "databases/garfield-data/annotation/"
  write.table(chr_df, file.path(output_dir, chr_string), sep = "", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
}

#Make a link file
link_file = data_frame(Index = c(0:7), annotation = names(all_qtl_pos), 
           Celltype = "Macrophage", "Tissue" = "Blood", 
           Type = c(rep("eQTL",4), rep("caQTL",4)), 
           Category = "QTL")
write.table(link_file, "databases/garfield-data/annotation/link_file.txt", quote = FALSE,
            sep = " ", row.names = FALSE, col.names = TRUE)




