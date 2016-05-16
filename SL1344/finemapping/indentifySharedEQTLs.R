library("GenomicRanges")
load_all("../seqUtils/")

#Import credible sets
credible_sets = readRDS("results/SL1344/eQTLs/rasqual_credible_sets.rds")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import RNA-Seq data
rna_list = readRDS("results/SL1344/combined_expression_data_covariates.rds")
min_pvalues_list = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df = plyr::ldply(min_pvalues_hits, .id = "condition_name") %>% 
  dplyr::group_by(gene_id) %>% dplyr::arrange(p_nominal)

#Create a df of SNP positions
snp_pos_df = vcf_file$snpspos %>% 
  dplyr::transmute(seqnames = chr, start = pos, end = pos, strand = "+", snp_id = snpid)

#Convert credible sets into gigantic data frame
credible_sets_df = purrr::map_df(credible_sets, ~purrr::map_df(., 
                                 ~dplyr::mutate(.,chr = as.character(chr))) %>% 
                                 dplyr::filter(chr != "X"), .id = "condition_name")

#Construct granges object
credible_sets_granges = credible_sets_df %>% 
  dplyr::transmute(condition_name, gene_id, snp_id, chr, seqnames = chr, start = pos, end = pos, strand = "+") %>% 
  dataFrameToGRanges()

#Find overlaps between credible sets
olaps = findOverlaps(credible_sets_granges)
queries = credible_sets_granges[queryHits(olaps),] %>% 
  elementMetadata() %>% as.data.frame() %>% 
  dplyr::transmute(master_condition = condition_name, master_id = gene_id, chr1 = chr) %>% 
  tbl_df()
subjects = credible_sets_granges[subjectHits(olaps),] %>% 
  elementMetadata() %>% as.data.frame() %>% 
  dplyr::transmute(dependent_condition = condition_name, dependent_id = gene_id, chr2 = chr) %>% 
  tbl_df()

#Identify all genes that share at least one SNP in their credible set
pairwise_shared = dplyr::bind_cols(queries, subjects) %>% 
  unique() %>% 
  dplyr::filter(chr1 == chr2) %>% 
  dplyr::filter(master_id != dependent_id) %>% 
  dplyr::select(master_id, master_condition, dependent_id, dependent_condition) %>% 
  unique()

#Calculate Jacccard scores for the credible sets
pairwise_shared_jaccard = purrr::by_row(pairwise_shared, credibleSetJaccard, credible_sets, .collate = "rows")
gene_pairs = dplyr::filter(pairwise_shared_jaccard, jaccard > 0.5) %>% 
  dplyr::group_by(master_id, dependent_id) %>%
  dplyr::arrange(p_nominal) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::group_by(master_id, dependent_id,snp_id) %>%
  dplyr::summarise(jaccard = max(jaccard)) %>%
  dplyr::ungroup()

#Convert shared peaks into clusters
clusters = constructClustersFromGenePairs(dplyr::select(gene_pairs, master_id, dependent_id))

#Use Coloc to test if the QTLs are indeed shared

#Test for GxGxE interactions
#Extract data
naive_ifng_rna = extractConditionFromExpressionList(c("naive","IFNg"), rna_list)
naive_sl1344_rna = extractConditionFromExpressionList(c("naive","SL1344"), rna_list)
naive_ifng_sl1344_rna = extractConditionFromExpressionList(c("naive","IFNg_SL1344"), rna_list)

#Specify models of interest
covariate_names = c("sex_binary", "ng_ul_mean","macrophage_diff_days","rna_auto", "max_purity_filtered", "harvest_stimulation_days",
                    "PEER_factor_1", "PEER_factor_2", "PEER_factor_3","PEER_factor_4", "PEER_factor_5","PEER_factor_6")

model0 = as.formula(paste("cqn ~ genotype + peak_type + condition_name + peak_type*condition_name + 
                    genotype*peak_type + genotype*condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
model1 = as.formula(paste("cqn ~genotype + peak_type + condition_name + peak_type*condition_name + 
                    genotype*peak_type + genotype*condition_name + genotype*condition_name*peak_type ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
peak_peak_interactions = purrr::by_row(gene_pairs, testThreewayInteraction, naive_ifng_rna$cqn, 
                                       naive_ifng_rna$sample_metadata, vcf_file, model0, model1, .collate = "rows", .to = "naive_vs_IFNg_pvalue")
naive_sl1344_ineractions = purrr::by_row(gene_pairs, testThreewayInteraction, naive_ifng_rna$cqn, 
                                       naive_ifng_rna$sample_metadata, vcf_file, model0, model1, .collate = "rows", .to = "naive_vs_SL1344_pvalue")
naive_ifng_sl1344_ineractions = purrr::by_row(gene_pairs, testThreewayInteraction, naive_ifng_rna$cqn, 
                                       naive_ifng_rna$sample_metadata, vcf_file, model0, model1, .collate = "rows", .to = "naive_vs_IFNg_SL1344_pvalue")

#Visualise a couple of examples
res = testThreewayInteraction(a[16,], naive_ifng_rna$cqn, naive_ifng_rna$sample_metadata, vcf_file, model0, model1, p_only = FALSE)
ggplot(res$data, aes(x = factor(genotype), y = cqn)) + geom_point() + geom_boxplot() + facet_grid(peak_type ~ condition_name)


