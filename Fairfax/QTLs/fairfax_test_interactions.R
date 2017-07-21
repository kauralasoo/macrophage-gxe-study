#Import SummarizedExperiment
se_fairfax = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")

#keep only sample that are present in all conditions
selected_samples = (colData(se_fairfax) %>% tbl_df2() %>% dplyr::filter(present_in_all))$sample_id
se_shared = se_fairfax[,selected_samples]

expression_mat = assays(se_shared)$exprs
sample_metadata = colData(se_shared) %>% tbl_df2() %>% 
  dplyr::rename(donor = donor_id)
gene_metadata = rowData(se_shared) %>% tbl_df2() %>% dplyr::mutate(gene_id = probe_id)

#Import QTLs
qtl_min_pvalues = readRDS("results/Fairfax/fairfax_qtl_min_pvalues.rds")

#Identify all qtl pairs
min_pvalue_hits = purrr::map(qtl_min_pvalues$shared, ~dplyr::filter(.,p_fdr < 0.1))
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(phenotype_id, p_nominal)
joint_pairs = dplyr::transmute(min_pvalues_df, gene_id = phenotype_id, snp_id, pheno_chr) %>% unique()

#Import genotypes
chr = "9"
vcf_file = seqUtils::gdsToMatrix(paste0("processed/Fairfax/geno_by_chr/",chr,".gds"))
chr_pairs = dplyr::filter(joint_pairs, pheno_chr == chr) %>%
  dplyr::select(gene_id, snp_id)

#Filter unique snps per probe
filtered_pairs = filterHitsR2(chr_pairs, vcf_file$genotypes, .8)

#Use a paired design to test for interaction
covariate_names = c("1")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name + (1|donor) ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype + (1|donor) ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(tbl_df(filtered_pairs), expression_mat, 
                                               sample_metadata, 
                                               vcf_file, formula_qtl, formula_interaction, 
                                               id_field_separator = "-", lme4 = TRUE)
interaction_df = postProcessInteractionPvalues(interaction_results, id_field_separator = "-")


