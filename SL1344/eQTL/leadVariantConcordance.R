library("devtools")
library("dplyr")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("../seqUtils/")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Calculate lead variant concordance between conditions
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
fastqtl_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")
fastqtl_100kb = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues_100kb.rds")


#Calculate lead variant sharing between pairs of conditions
concordance_matrix = calculatePairwiseConcordance(rasqual_min_pvalues, 1e-5, vcf_file$genotypes)
fastqtl_concordance_matrix = calculatePairwiseConcordance(fastqtl_min_pvalues, ~p_fdr < 0.01, vcf_file$genotypes)
fastqtl_100kb_concordance_matrix = calculatePairwiseConcordance(fastqtl_100kb, ~p_fdr < 0.01, vcf_file$genotypes)

#Rasqual vs fastQTL concordance
calculatePairwiseConcordance(list(rasqual = rasqual_min_pvalues$naive, fastqtl = fastqtl_min_pvalues$naive), 
                             1e-5, vcf_file$genotypes)


#Look at splicing QTLs
salmon_min_pvalues = readRDS("results/SL1344/salmon/salmon_qtl_hits.rds")
salmon_event_concordance = calculatePairwiseConcordance(salmon_min_pvalues, ~p_fdr < 0.01, vcf_file$genotypes)

#Aggregare reviseAnnotations at gene level
salmon_gene_pvalues = purrr::map(salmon_min_pvalues, ~dplyr::arrange(., ensembl_gene_id, p_nominal) %>% 
             dplyr::group_by(ensembl_gene_id) %>% 
             dplyr::filter(row_number() == 1) %>%
             dplyr::ungroup() %>%
             dplyr::mutate(gene_id = ensembl_gene_id))
salmon_gene_concordance = calculatePairwiseConcordance(salmon_gene_pvalues, ~p_fdr < 0.01, vcf_file$genotypes)

salmon_ensembl_pvalues = readRDS("results/SL1344/salmon/salmon_ensembl_qtl_hits.rds")
salmon_ensembl_concordance = calculatePairwiseConcordance(salmon_ensembl_pvalues, ~p_fdr < 0.01, vcf_file$genotypes)

calculatePairwiseConcordance(list(rasqual = salmon_gene_pvalues$naive, fastqtl = salmon_ensembl_pvalues$naive), 
                             filter_formula = ~p_fdr < 0.01, vcf_file$genotypes)
                                        

#Concordane betweeb tQTLs and eQTLs
calculatePairwiseConcordance(list(rasqual = fastqtl_100kb$naive, fastqtl = salmon_gene_pvalues$naive), 
                            filter_formula = ~p_fdr < 0.01, vcf_file$genotypes)






