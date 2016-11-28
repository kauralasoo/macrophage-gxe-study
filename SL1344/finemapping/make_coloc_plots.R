library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("ggplot2")
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

#Make colocalisation plots
qtl_df = dplyr::filter(rasqual_min_pvalues$naive, gene_id == "ENSG00000179344") %>% 
  dplyr::select(gene_id, snp_id) 

results = colocMolecularQTLs(qtl_df ,qtl_summary_path = qtlResults()$rna_rasqual$naive, 
                             gwas_summary_path = "databases/GWAS/GWAS/Alzheimers_disease_Lambert_2013_NatGen_GWAS_meta_stage1.sorted.txt.gz"
                             ,GRCh37_variants, GRCh38_variants, qtl_type = "rasqual")

#Join data together
trait_df = purrr::map_df(results$data, ~dplyr::select(.,chr, pos, p_nominal), .id = "trait") %>%
  dplyr::mutate(log10p = -log(p_nominal, 10))

#Make a plot
ggplot(trait_df, aes(x = pos, y = log10p)) + 
  geom_point() + facet_wrap(~trait, ncol = 1, scales = "free_y")

