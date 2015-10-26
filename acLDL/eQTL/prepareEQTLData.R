library("devtools")
library("dplyr")
load_all("../seqUtils/")
library("SNPRelate")

#Import data
expression_dataset = readRDS("results/acLDL/acLDL_combined_expression_data.rds") #expression data
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
acLDL_meta = read.table("macrophage-gxe-study/data/sample_lists/acLDL_line_data.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(line_id, replicate) %>%
  tidyr::separate(line_id, c("donor"), remove = FALSE, extra ="drop") %>%
  unique()
gene_id_name_map = dplyr::select(expression_dataset$gene_metadata, gene_id, gene_name)

#Filter out some samples and discard replicates
design = dplyr::filter(expression_dataset$design, !(donor == "mijn")) %>% tbl_df() %>% #Should not have received from cGAP
  dplyr::filter(!(donor == "xegx")) %>% #Very strong outlier on the heatmap
  dplyr::left_join(acLDL_meta, by = "donor")
sample_meta = dplyr::left_join(design, line_metadata, by = c("donor", "replicate", "line_id")) %>%
  dplyr::mutate(condition_name = factor(condition, levels = c("AcLDL","Ctrl")))

#Construct separate design matrices
cond_A_design = dplyr::filter(sample_meta, condition == "Ctrl")
cond_B_design = dplyr::filter(sample_meta, condition == "AcLDL")
design_list = list(Ctrl = cond_A_design, AcLDL = cond_B_design)

#### EXPRESSION ####
#Filter expression data by samples
exprs_cqn = expression_dataset$exprs_cqn[,design$sample_id]
exprs_counts = expression_dataset$exprs_counts[,design$sample_id]

#Define expressed genes
mean_expression = calculateMean(exprs_cqn, as.data.frame(sample_meta), "condition_name")
expressed_genes = names(which(apply(mean_expression, 1, max) > 0))

#Set up the genepos data.frame
genepos = dplyr::filter(expression_dataset$gene_metadata, gene_id %in% expressed_genes) %>% 
  dplyr::transmute(chr = chromosome_name, left = start_position, right = end_position, geneid = gene_id, score = 1000, strand) %>%
  dplyr::mutate(strand = ifelse(strand == 1, "+","-")) %>%
  as.data.frame() %>%
  dplyr::filter( !(chr %in% c("MT","Y")) ) %>% #Remove genes on MT and Y chromosomes
  dplyr::arrange(chr, left, right)

#Filter expression data by min expression
exprs_cqn = expression_dataset$exprs_cqn[genepos$geneid, design$sample_id]
exprs_counts = expression_dataset$exprs_counts[genepos$geneid, design$sample_id]

#Set up expression data for each condition
condA_exp = extractSubset(dplyr::filter(sample_meta, condition == "Ctrl"), exprs_cqn)
condB_exp = extractSubset(dplyr::filter(sample_meta, condition == "AcLDL"), exprs_cqn)
exprs_cqn_list = list(Ctrl = condA_exp, AcLDL = condB_exp)

#### COVARIATES ####
#Construct covariate matrix
covariates = dplyr::mutate(sample_meta, sex = ifelse(gender == "male",1,0)) %>%
  dplyr::select(sample_id, sex)

#Save expression data for PEER and run PEER outside of R
savePEERData(exprs_cqn_list, "results/acLDL/PEER/input/")
#Run PEER .....

#Import PEER results back into R
cond_A_factors = importPEERFactors("results/acLDL/PEER/Ctrl_10/factors.txt", cond_A_design)
cond_B_factors = importPEERFactors("results/acLDL/PEER/AcLDL_10/factors.txt", cond_B_design)
peer_factors = rbind(cond_A_factors,cond_B_factors) %>%
  dplyr::semi_join(design, by = "sample_id")

#Merge PEER facotrs with covariates
peer_covariates = dplyr::left_join(covariates, peer_factors, by = "sample_id") %>% as.data.frame()
rownames(peer_covariates) = peer_covariates$sample_id
peer_covariates = t(peer_covariates[,-1])

#Set up covatiates for each conditon
condA_covs = extractSubset(dplyr::filter(sample_meta, condition == "Ctrl"), peer_covariates)
condB_covs = extractSubset(dplyr::filter(sample_meta, condition == "AcLDL"), peer_covariates)
covariates_list = list(Ctrl = condA_covs, AcLDL = condB_covs)

#### Dataset object ####
#Construct a single dataset object that can be reused by many different analysis
eqtl_dataset = list(
  exprs_cqn_list = exprs_cqn_list,
  covariates_list = covariates_list,
  genepos = genepos,
  gene_metadata = expression_dataset$gene_metadata,
  sample_metadata = sample_meta,
  covariates = t(peer_covariates),
  exprs_cqn = exprs_cqn
)
saveRDS(eqtl_dataset, "results/acLDL/acLDL_eqtl_data_list.rds")


### FastQTL ####
#Export expression data for fastQTL analysis
eqtl_dataset = readRDS("results/acLDL/acLDL_eqtl_data_list.rds")

#Extract data from list
donor_genotype_map = dplyr::filter(eqtl_dataset$sample_metadata, condition == "Ctrl") %>% dplyr::select(donor, genotype_id)
genepos = dplyr::select(eqtl_dataset$genepos, chr, left, right, geneid)
colnames(genepos)[1] = "#chr"

#Make an updated list of covariates
fastqtl_covariates_list = lapply(eqtl_dataset$covariates_list, function(cov_mat, donor_genotype_map){
  res = extractSubset(donor_genotype_map, as.data.frame(cov_mat), "donor", "genotype_id") %>% 
    dplyr::mutate(id = rownames(cov_mat)) %>% 
    dplyr::select(id, everything())
  return(res)
}, donor_genotype_map) %>%
  lapply(., function(x){x[1:7,]}) #Keep the first 7 covariates

#Make an updated list of gene expression
fastqtl_expression_list = lapply(eqtl_dataset$exprs_cqn_list, function(exp_mat, donor_genotype_map, genepos){
  res = extractSubset(donor_genotype_map, as.data.frame(exp_mat), "donor","genotype_id") %>%
    dplyr::mutate(geneid = rownames(exp_mat)) %>%
    dplyr::select(geneid, everything()) %>%
    dplyr::left_join(genepos, ., by = "geneid")
  return(res)
}, donor_genotype_map, genepos)

#Save matrices to disk
saveFastqtlMatrices(fastqtl_covariates_list, "results/acLDL/fastqtl/input/", file_suffix = "covariates")
saveFastqtlMatrices(fastqtl_expression_list, "results/acLDL/fastqtl/input/", file_suffix = "expression")
write.table(donor_genotype_map$genotype_id, "results/acLDL/fastqtl/input/genotype_list.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Construct chunks table
chunks_matrix = data.frame(chunk = seq(1:200), n = 200)
write.table(chunks_matrix, "results/acLDL/fastqtl/input/chunk_table.txt", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")



