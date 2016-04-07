library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")


#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Load p-values from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_min_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df = plyr::ldply(rasqual_min_hits, .id = "condition_name") %>% dplyr::arrange(p_nominal)
joint_pairs = dplyr::select(min_pvalue_df, gene_id, snp_id) %>% unique() 

#List of RNA-seq tabix files
rna_list = list(naive = "databases/SL1344/naive_500kb.sorted.txt.gz",
                IFNg = "databases/SL1344/IFNg_500kb.sorted.txt.gz",
                SL1344 = "databases/SL1344/SL1344_500kb.sorted.txt.gz",
                IFNg_SL1344 = "databases/SL1344/IFNg_SL1344_500kb.sorted.txt.gz")
                

#List of ATAC tabix files
atac_tabix_list = list(naive = "../macrophage-chromatin/results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz",
                       IFNg = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz",
                       SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz",
                       IFNg_SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz")

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Fetch top genes
naive_list = idVectorToList(rasqual_min_hits$naive$gene_id)
naive_cs = purrr::map(naive_list, ~tabixFetchGenesQuick(.,rna_list$naive, combined_expression_data$gene_metadata, cis_window = 5e5)[[1]] %>% 
                        addAssociationPosterior(69) %>% constructCredibleSet(0.999))

#Find hits for naive QTLs
naive_top_hits = purrr::map(naive_list, ~tabixFetchGenesQuick(.,rna_list$naive, combined_expression_data$gene_metadata, cis_window = 5e5)[[1]] %>% 
                        dplyr::arrange(-chisq) %>% head(50))
naive_annotated = purrr::map(naive_top_hits, ~annotateCredibleSet(.,atac_list$gene_metadata, atac_tabix_list$naive))
naive_filtered = purrr::map(naive_annotated, ~dplyr::filter(., overlap_peak_id == assoc_peak_id))


#Overlap
naive_credible_set = constructCredibleSet(naive_pvalues, threshold = 0.99)
naive_cs_annotated = annotateCredibleSet(naive_credible_set, atac_list$gene_metadata, atac_tabix_list$naive)

ifng_sl1344_cs = constructCredibleSet(IFNg_SL1344_pvalues, threshold = 0.99)
ifng_sl1344_cs_annotated = annotateCredibleSet(ifng_sl1344_cs, atac_list$gene_metadata, atac_tabix_list$IFNg_SL1344)

