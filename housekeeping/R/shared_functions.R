#Helper functions
enrichFastQTLPvalues <- function(fastqtl_pvalues, gene_id_name_map){
  res = tbl_df(fastqtl_pvalues) %>%
    dplyr::arrange(p_beta) %>%
    dplyr::filter(!is.na(p_beta)) %>%
    dplyr::mutate(p_beta_log10 =-log(p_beta,10)) %>% 
    dplyr::mutate(p_expected = -log(c(1:length(p_beta))/length(p_beta),10)) %>% 
    dplyr::left_join(gene_id_name_map, by = "gene_id") %>%
    dplyr::select(gene_id, gene_name, everything()) %>%
    dplyr::mutate(qvalue = qvalue::qvalue(p_beta)$qvalues)
  return(res)
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