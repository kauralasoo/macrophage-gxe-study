library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-chromatin/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")
library("purrr")

#### Import data ####
#Load the raw eQTL dataset
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(atac_list$gene_metadata, gene_id, gene_name)

#Load min p-values with optimized genotypes
#Load p-values from disk
rasqual_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
rasqual_min_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df = plyr::ldply(rasqual_min_hits, .id = "condition_name") %>% dplyr::arrange(p_nominal)

#Min p-values with fixed genotypes
#Load p-values from disk
rasqual_min_pvalues_fixed = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.fixed_gt.rds")
rasqual_min_hits_fixed = lapply(rasqual_min_pvalues_fixed, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df_fixed = plyr::ldply(rasqual_min_hits_fixed, .id = "condition_name") %>% dplyr::arrange(p_nominal)

#List of RNA-seq tabix files
atac_tbx_list = list(naive = "results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz",
                IFNg = "results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz",
                SL1344 = "results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz",
                IFNg_SL1344 = "results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz")
atac_tbx_list_fixed = list(naive = "results/ATAC/rasqual/output_fixed_gt/naive_100kb/naive_100kb.sorted.txt.gz",
                IFNg = "results/ATAC/rasqual/output_fixed_gt/IFNg_100kb/IFNg_100kb.sorted.txt.gz",
                SL1344 = "results/ATAC/rasqual/output_fixed_gt/SL1344_100kb/SL1344_100kb.sorted.txt.gz",
                IFNg_SL1344 = "results/ATAC/rasqual/output_fixed_gt/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz")

#Import the VCF file
vcf_file = readRDS("../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

naive_list = idVectorToList(rasqual_min_hits$naive$gene_id)[1:10]
naive_cs = purrr::map(naive_list, ~tabixFetchGenesQuick(.,atac_tbx_list$naive, atac_list$gene_metadata, cis_window = 5e4)[[1]] %>%
                        dplyr::arrange(p_nominal) %>% addR2FromLead(vcf_file$genotypes) %>% dplyr::filter(R2 > 0.8))

naive_cs_fixed = purrr::map(naive_list, ~tabixFetchGenesQuick(.,atac_tbx_list_fixed$naive, atac_list$gene_metadata, cis_window = 5e4)[[1]] %>%
                        dplyr::arrange(p_nominal) %>% addR2FromLead(vcf_file$genotypes) %>% dplyr::filter(R2 > 0.8))

ranges = constructGeneRanges(rasqual_min_hits$naive[1:100,], atac_list$gene_metadata, 5e4)
res = tabixFetchGenes(ranges, atac_tbx_list$naive)
cs_list = purrr::map(res, ~dplyr::arrange(.,p_nominal) %>% addR2FromLead(vcf_file$genotypes) %>% dplyr::filter(R2 > 0.8))

#Fixed gt
res_fixed = tabixFetchGenes(ranges, atac_tbx_list_fixed$naive)
cs_list_fixed = purrr::map(res_fixed, ~dplyr::arrange(.,p_nominal) %>% addR2FromLead(vcf_file$genotypes) %>% dplyr::filter(R2 > 0.8))


cs_list_chisq = purrr::map(res, ~dplyr::filter(., chisq > max(chisq)*0.9))


