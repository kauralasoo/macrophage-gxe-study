library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import eQTLs
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")

#Import GRCh37 variant information
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")

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






