constructDesignMatrix_acLDL <- function(sample_ids){
  #Construct a design matrix from a list of sample ids
  
  design = data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  design = design %>% 
    tidyr::separate(sample_id, into = c("donor", "time_point", "condition"), remove = FALSE, extra = "merge") %>% 
    dplyr::mutate(donor = tolower(donor)) %>%
    dplyr::mutate(condition_name = condition)
  rownames(design)= design$sample_id
  return(design)
}