library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import the VCF file
vcf_file = readRDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.rds")

#Import GWAS olaps
gwas_olaps = readRDS("acLDL_figures/tables/GWAS_coloc_hits.rds")

#Import interactions test results
trQTL_interactions = readRDS("results/acLDL/trQTLs/trQTL_interaction_results.rds")

#Import eQTL interaction test results
#NOTE: there seems to be a problem with a difference between RASQUAL and LM eQTLs
interaction_df = readRDS("results/acLDL/eQTLs/lme4_qtltools_interaction_results.rds")
interaction_hits = dplyr::filter(interaction_df, p_fdr < 0.01) %>%
  dplyr::transmute(phenotype_id = gene_id, qtl_snp_id = snp_id)
interaction_eqtl_olaps = dplyr::left_join(interaction_hits, gwas_olaps$featureCounts, by = "phenotype_id") %>%
  dplyr::filter(!is.na(trait)) %>%
  dplyr::group_by(phenotype_id, snp_id, condition_name, qtl_snp_id, trait) %>%
  dplyr::mutate(R2 = calculatePairR2(qtl_snp_id, snp_id, vcf_file$genotypes)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(R2 > 0.8)
write.table(interaction_eqtl_olaps, "acLDL_figures/tables/GWAS_coloc_eQTL_interactions.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)


#trQTL interaction hits
trQTL_hits = purrr::map(trQTL_interactions, ~dplyr::filter(.,p_fdr < 0.01) %>% 
             dplyr::transmute(phenotype_id = gene_id, qtl_snp_id = snp_id))

trQTL_iteraction_olaps = purrr::map2(trQTL_hits, gwas_olaps[1:3], ~dplyr::left_join(.x, .y, by = "phenotype_id") %>%
              dplyr::filter(!is.na(trait))) %>%
  purrr::map_df(identity, .id = "annotation") %>%
  dplyr::group_by(phenotype_id, snp_id, condition_name, qtl_snp_id, trait) %>%
  dplyr::mutate(R2 = calculatePairR2(qtl_snp_id, snp_id, vcf_file$genotypes)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(R2 > 0.8)
write.table(trQTL_iteraction_olaps, "acLDL_figures/tables/GWAS_coloc_trQTL_interactions.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)





