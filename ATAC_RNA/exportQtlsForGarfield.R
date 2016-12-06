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
                                                     GRCh37_variants, atac_data$gene_metadata, coords, cis_dist = 1e5)
  
  #Export annotations
  output_dir = "databases/garfield-data/annotation_caqtl/"
  write.table(garfield_results, file.path(output_dir, chr_string), sep = "", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  
}

#Make a link file
link_file = data_frame(Index = c(0:3), annotation = names(rasqual_min_pvalues), 
                       Celltype = "Macrophage", "Tissue" = "Blood", 
                       Type = rep("caQTL",4), 
                       Category = "QTL")
write.table(link_file, "databases/garfield-data/annotation/link_file.txt", quote = FALSE,
            sep = " ", row.names = FALSE, col.names = TRUE)
