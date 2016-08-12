#This script performs the timecomsuming bit of fitting mixed-effects
#models on normalized gene expression data
library("dplyr")
library("lme4")
library("devtools")
library("purrr")
load_all("../seqUtils/")

batch_id = NULL

####### Get batch id from disk ######
f <- file("stdin")
open(f)
batch_id = readLines(f) %>% as.numeric()
close(f)
####### END #######

#consturct VarComp metadata
varCompMeta <- function(sample_metadata){
  metadata = sample_metadata %>% 
    dplyr::select(sample_id, line_id, SL1344, IFNg, rna_auto, ng_ul_mean, 
                  salmonella_date, chemistry, sex, mean_purity_filtered, passage_diff,
                  macrophage_diff_days, library_pool, rna_extraction) %>%
    dplyr::mutate(chemistry = ifelse(chemistry == "V3", "V3", "V4")) %>%
    dplyr::mutate(purity_bins = cut(mean_purity_filtered, breaks = c(0,0.95,0.975,1))) %>%
    dplyr::mutate(concentration_bins = cut(ng_ul_mean, breaks = c(0,100,200,300,500))) %>%
    dplyr::mutate(passage_diff_bins = cut(passage_diff, breaks = c(0,25,35,45,60))) %>%
    dplyr::mutate(macrophage_diff_bins = cut(macrophage_diff_days, breaks = c(0,30,40,50,60,100)))
  return(metadata)
}

#Import eQTL data list
data_list = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#Extract separate lists for each condition
condition_names = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
data_by_conditons = lapply(condition_names, extractConditionFromExpressionList, data_list)

#Extract expression matrix and metadata matrices from condtion list
cqn_conditions = purrr::map(data_by_conditons, ~.$cqn)
metadata_conditions = purrr::map(data_by_conditons, ~.$sample_metadata)
varcomp_meta_conditons = purrr::map(metadata_conditions, varCompMeta)

#### Split genes into batches ####
gene_ids = rownames(data_list$cqn)
batches = splitIntoBatches(length(gene_ids), 50)
if(!is.null(batch_id)){
  selection = batches == batch_id
  gene_ids = gene_ids[selection]
}

#Prepare data for lmer analysis
gene_id_list = idVectorToList(gene_ids)

#Make gene lists for each condition
gene_data_list = purrr::map2(cqn_conditions, varcomp_meta_conditons, function(x,y,gene_id_list){
  purrr::map(gene_id_list, ~constructGeneData(., x, y))
}, gene_id_list)

#Fit the extended model
model_extended <- function(model_data){
  mod = lme4::lmer(value ~ (1|rna_auto) + (1|salmonella_date) + (1|sex) + (1|chemistry) + 
                     (1|purity_bins) + (1|concentration_bins) + (1|passage_diff_bins) + 
                     (1|macrophage_diff_bins) + (1|library_pool) +
                     (1|rna_extraction), model_data)
  return(mod) 
}

#Calculate variance explained for each gene in each condition
variance_list_list = purrr::map(gene_data_list, ~purrr::map(., estimateVarianceExplained, model_extended))
var_table = purrr::map_df(variance_list_list, ~purrr::map_df(., identity, .id = "gene_id"), .id = "condition_name")

#Save output from each batch
if(!is.null(batch_id)){
  output_file = file.path("results/SL1344/varComp/", paste0("varCompByCondition_batch_",batch_id, ".txt"))
  write.table(var_table, output_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

