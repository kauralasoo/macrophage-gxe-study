library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("devtools")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import eQTLs
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
fastqtl_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")

#Import old and new variant coordinates
GRCh38_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")

#Filter gene
qtl_df = dplyr::filter(rasqual_min_pvalues$naive, p_fdr < 0.1) %>% head(10) %>% dplyr::select(gene_id, snp_id)

#Test for coloc
coloc_res = purrr::by_row(qtl_df, ~colocMolecularQTLs(.,qtl_summary_path = qtlResults()$rna_rasqual$naive, 
    gwas_summary_path = "databases/GWAS/GWAS/Alzheimers_disease_Lambert_2013_NatGen_GWAS_meta_stage1.sorted.txt.gz"
    ,GRCh37_variants, GRCh38_variants, qtl_type = "rasqual")$summary, .collate = "rows")


    
  

