#Construct design matrix
constructDesignMatrix_ATAC <- function(sample_names){
  #Prespecify conditon names
  condition_names = data.frame(condition_char = c("A","B","C","D"), 
                               condition_name = c("naive", "IFNg", "SL1344", "IFNg_SL1344"), 
                               stringsAsFactors = FALSE)
  
  design = data_frame(sample_id = sample_names) %>% 
    tidyr::separate(sample_id, c("donor", "condition_char", "assay"), sep = "_",remove = FALSE) %>%
    dplyr::left_join(condition_names, by = "condition_char") %>%
    dplyr::select(sample_id, donor, condition_char, condition_name) %>%
    dplyr::mutate(SL1344 = ifelse(condition_char %in% c("C","D"), "infected", "control")) %>%
    dplyr::mutate(IFNg = ifelse(condition_char %in% c("B","D"), "primed", "naive"))
  return(design)
}

constructDesignMatrix_Ivashkiv <- function(sample_names){
  #Prespecify conditon names
  condition_names = data.frame(condition_char = c("A","B","C","D"), 
                               condition_name = c("naive", "IFNg", "LPS", "IFNg_LPS"), 
                               stringsAsFactors = FALSE)
  
  design = data_frame(sample_id = sample_names) %>% 
    tidyr::separate(sample_id, c("antibody", "replicate", "condition_char"), sep = "_",remove = FALSE) %>%
    dplyr::left_join(condition_names, by = "condition_char") %>%
    dplyr::mutate(LPS = ifelse(condition_char %in% c("C","D"), "stimulated", "control")) %>%
    dplyr::mutate(IFNg = ifelse(condition_char %in% c("B","D"), "primed", "naive"))
  return(design)
}