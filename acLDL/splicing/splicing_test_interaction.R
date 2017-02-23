library("devtools")
library("dplyr")
library("lme4")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Explore the principal components of splicing data

#Import datasets
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
se_revised = readRDS("results/acLDL/acLDL_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/acLDL/acLDL_leafcutter_counts.rds")

#Calculate normalised abundances
ensembl_abundances = assays(se_ensembl)$tpm_ratios %>% 
  replaceNAsWithRowMeans() %>% 
  quantileNormaliseRows()
revised_abundances = assays(se_revised)$tpm_ratios %>% 
  replaceNAsWithRowMeans() %>% 
  quantileNormaliseRows()
leafcutter_abundances = assays(se_leafcutter)$tpm_ratios %>% 
  replaceNAsWithRowMeans() %>% 
  quantileNormaliseRows()

#Peform PCA
ensembl_pca = performPCA(ensembl_abundances,  tbl_df2(colData(se_ensembl)), 
                         n_pcs = 10, feature_id = "sample_id")$pca_matrix
revised_pca = performPCA(revised_abundances,  tbl_df2(colData(se_revised)), 
                         n_pcs = 10, feature_id = "sample_id")$pca_matrix 
leafcutter_pca = performPCA(leafcutter_abundances,  tbl_df2(colData(se_leafcutter)), 
                         n_pcs = 10, feature_id = "sample_id")$pca_matrix 

#Make some plots
ggplot(ensembl_pca, aes(x = PC1, y = PC2, color = factor(rna_submit), shape = condition_name)) + geom_point()
ggplot(revised_pca, aes(x = PC1, y = PC2, color = factor(condition_name), shape = condition_name)) + geom_point()
ggplot(leafcutter_pca, aes(x = PC1, y = PC2, color = factor(condition_name), shape = condition_name)) + geom_point()

#Extract covariates
ensembl_covariates = c("PC1", "PC2", "PC3", "PC6")
revised_covariates = c("PC1", "PC2", "PC3", "PC4")
leafcutter_covariates = c("PC1", "PC3", "PC4", "PC5", "PC6")

#Import the VCF file
vcf_file = readRDS("genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.rds")


#### Ensembl interactions ####
#Import p-values
ensembl_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/AcLDL.permuted.txt.gz"))

#Identify mininal p-values
min_pvalue_hits = purrr::map(ensembl_pvalues, ~dplyr::filter(.,p_fdr < 0.1))
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(phenotype_id, p_nominal) %>%
  dplyr::mutate(gene_id = phenotype_id)

#Filter gene-SNP pairs by R2
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

##### Ctrl vs AcLDL (paired design with lme4) #####
covariate_names = ensembl_covariates
formula_qtl = as.formula(paste("expression ~ genotype + condition_name + (1|donor) ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype + (1|donor) ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(filtered_pairs, ensembl_abundances, ensembl_pca, 
                                               vcf_file, formula_qtl, formula_interaction, id_field_separator = "-", lme4 = TRUE)
ensembl_df = postProcessInteractionPvalues(interaction_results, id_field_separator = "-")


#### Revised interactions ####
#Import p-values
revised_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/reviseAnnotations/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/fastqtl_output/reviseAnnotations/AcLDL.permuted.txt.gz"))

#Identify mininal p-values
min_pvalue_hits = purrr::map(revised_pvalues, ~dplyr::filter(.,p_fdr < 0.1))
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(phenotype_id, p_nominal) %>%
  dplyr::mutate(gene_id = phenotype_id)

#Filter gene-SNP pairs by R2
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

##### Ctrl vs AcLDL (paired design with lme4) #####
covariate_names = revised_covariates
formula_qtl = as.formula(paste("expression ~ genotype + condition_name + (1|donor) ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype + (1|donor) ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(filtered_pairs, revised_abundances, revised_pca, 
                                               vcf_file, formula_qtl, formula_interaction, id_field_separator = "-", lme4 = TRUE)
revised_df = postProcessInteractionPvalues(interaction_results, id_field_separator = "-")


#### LeafCutter interactions ####
leafcutter_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/leafcutter/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/fastqtl_output/leafcutter/AcLDL.permuted.txt.gz"))

#Identify mininal p-values
min_pvalue_hits = purrr::map(leafcutter_pvalues, ~dplyr::filter(.,p_fdr < 0.1))
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(phenotype_id, p_nominal) %>%
  dplyr::mutate(gene_id = phenotype_id)

#Filter gene-SNP pairs by R2
joint_pairs = dplyr::select(min_pvalues_df, gene_id, snp_id) %>% unique()
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

##### Ctrl vs AcLDL (paired design with lme4) #####
covariate_names = leafcutter_covariates
formula_qtl = as.formula(paste("expression ~ genotype + condition_name + (1|donor) ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype + (1|donor) ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(filtered_pairs, leafcutter_abundances, leafcutter_pca, 
                                               vcf_file, formula_qtl, formula_interaction, id_field_separator = "-", lme4 = TRUE)
leafcutter_df = postProcessInteractionPvalues(interaction_results, id_field_separator = "-")


