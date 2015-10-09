constructDesignMatrix_SL1344 <- function(sample_ids){
  #Construct a design matrix from a list of sample ids
  
  design = data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  condition_names = data.frame(condition = c("A","B","C","D"), 
                               condition_name = c("naive", "IFNg", "SL1344", "IFNg_SL1344"), 
                               stringsAsFactors = FALSE)
  design = design %>% 
    tidyr::separate(sample_id, into = c("donor", "condition", "replicate"), extra = "drop", remove = FALSE) %>%
    dplyr::mutate(replicate = ifelse(is.na(replicate), 1, replicate)) %>% #if no replicate in file name then set it to 1
    dplyr::mutate(replicate = as.numeric(replicate)) %>%
    dplyr::left_join(condition_names, by = "condition") %>%
    dplyr::mutate(SL1344 = ifelse(condition %in% c("C","D"), "infected", "control")) %>%
    dplyr::mutate(IFNg = ifelse(condition %in% c("B","D"), "primed", "naive"))
  rownames(design)= design$sample_id
  return(design)
}

testInteraction <- function(gene_id, snp_id, eqtl_data_list){
  #Test for interaction
  sample_meta = dplyr::select(eqtl_data_list$sample_metadata, sample_id, donor, condition_name)
  cov_mat = dplyr::mutate(as.data.frame(eqtl_data_list$covariates), sample_id = rownames(eqtl_data_list$covariates))
  
  #Extract data
  exp_data = data_frame(sample_id = colnames(eqtl_data_list$exprs_cqn), expression = eqtl_data_list$exprs_cqn[gene_id,])
  geno_data = data_frame(donor = colnames(eqtl_data_list$genotypes), genotype = eqtl_data_list$genotypes[snp_id,])
  
  sample_data = dplyr::left_join(sample_meta, cov_mat, by = "sample_id") %>%
    dplyr::left_join(exp_data, by = "sample_id") %>%
    dplyr::left_join(geno_data, by = "donor")
  
  no_interaction = lm(expression~genotype + condition_name + sex + ng_ul_mean + diff_days + PEER_factor_1 + PEER_factor_2 + PEER_factor_3 + PEER_factor_4, as.data.frame(sample_data))
  interaction = lm(expression~genotype + condition_name + condition_name:genotype + sex + ng_ul_mean + diff_days + PEER_factor_1 + PEER_factor_2 + PEER_factor_3 + PEER_factor_4, as.data.frame(sample_data))
  res = anova(no_interaction, interaction)
  return(list(anova = res, interaction_model = interaction))
}

testMultipleInteractions <- function(snps_df, eqtl_data_list){
  #Plot eQTL results for a list of gene and SNP pairs.
  result = list()
  for(i in 1:nrow(snps_df)){
    gene_id = snps_df[i,]$gene_id
    snp_id = snps_df[i,]$snp_id
    print(gene_id)
    test = testInteraction(gene_id, snp_id, eqtl_data_list)
    result[[gene_id]] = test
  }
  return(result)
}

constructGeneData <- function(gene_id, snp_id, eqtl_data_list){
  #Test for interaction
  sample_meta = dplyr::select(eqtl_data_list$sample_metadata, sample_id, donor, condition_name)
  cov_mat = dplyr::mutate(as.data.frame(eqtl_data_list$covariates), sample_id = rownames(eqtl_data_list$covariates))
  
  #Extract data
  exp_data = data_frame(sample_id = colnames(eqtl_data_list$exprs_cqn), expression = eqtl_data_list$exprs_cqn[gene_id,])
  geno_data = data_frame(donor = colnames(eqtl_data_list$genotypes), genotype = eqtl_data_list$genotypes[snp_id,])
  
  sample_data = dplyr::left_join(sample_meta, cov_mat, by = "sample_id") %>%
    dplyr::left_join(exp_data, by = "sample_id") %>%
    dplyr::left_join(geno_data, by = "donor")
  
  return(sample_data)
}
