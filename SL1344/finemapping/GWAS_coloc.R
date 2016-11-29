library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("devtools")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import old and new variant coordinates
GRCh38_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")

#Import gene names
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

####### IBD ######
#Import eQTLs
qtl_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")

#Prefilter potential GWAS overlaps
gwas_prefix = "databases/GWAS/GWAS/Inflammatory_bowel_disease_Liu_2015_NatGen_Immunochip"
qtl_df_list = prefilterColocCandidates(qtl_min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.1, overlap_dist = 1e5, gwas_thresh = 1e-5)

#Test for coloc
coloc_res_list = purrr::map2(qtl_df_list, qtlResults()$rna_fastqtl, ~colocMolecularQTLsByRow(.x, qtl_summary_path = .y, 
                   gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants = GRCh37_variants, 
                   GRCh38_variants = GRCh38_variants, qtl_type = "fastqtl", N_qtl = 84, cis_dist = 1e5))

#Identify all unique coloc hits
coloc_hits = purrr::map_df(coloc_res_list, ~dplyr::filter(.,PP.H4.abf > 0.8) %>% 
                                   dplyr::left_join(gene_name_map, by = "gene_id"), .id = "condition_name")
write.table(coloc_hits, "results/SL1344/coloc/coloc_lists/IBD.coloc.txt", sep = "\t", quote = FALSE, row.names = FALSE)

####### SLE ######
#Import eQTLs
qtl_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")

#Prefilter potential GWAS overlaps
gwas_prefix = "databases/GWAS/GWAS/Systemic_lupus_erythematosus_Bentham_2015_NatGen_GWAS"
qtl_df_list = prefilterColocCandidates(qtl_min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.1, overlap_dist = 1e5, gwas_thresh = 1e-5)

#Test for coloc
coloc_res_list = purrr::map2(qtl_df_list, qtlResults()$rna_fastqtl, ~colocMolecularQTLsByRow(.x, qtl_summary_path = .y, 
            gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants = GRCh37_variants, 
            GRCh38_variants = GRCh38_variants, qtl_type = "fastqtl", N_qtl = 84, cis_dist = 1e5))

#Identify all unique coloc hits
coloc_hits = purrr::map_df(coloc_res_list, ~dplyr::filter(.,PP.H4.abf > 0.8) %>% 
                             dplyr::left_join(gene_name_map, by = "gene_id"), .id = "condition_name")
write.table(coloc_hits, "results/SL1344/coloc/coloc_lists/SLE.coloc.txt", sep = "\t", quote = FALSE, row.names = FALSE)

###### AZ ######
#Import eQTLs
qtl_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")

#Prefilter potential GWAS overlaps
gwas_prefix = "databases/GWAS/GWAS/Alzheimers_disease_Lambert_2013_NatGen_GWAS_meta_stage1"
qtl_df_list = prefilterColocCandidates(qtl_min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.1, overlap_dist = 1e5, gwas_thresh = 1e-5)

#Test for coloc
coloc_res_list = purrr::map2(qtl_df_list, qtlResults()$rna_fastqtl, ~colocMolecularQTLsByRow(.x, qtl_summary_path = .y, 
                                gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants = GRCh37_variants, 
                                GRCh38_variants = GRCh38_variants, qtl_type = "fastqtl", N_qtl = 84, cis_dist = 1e5))

#Identify all unique coloc hits
coloc_hits = purrr::map_df(coloc_res_list, ~dplyr::filter(.,PP.H4.abf > 0.8) %>% 
                             dplyr::left_join(gene_name_map, by = "gene_id"), .id = "condition_name")
write.table(coloc_hits, "results/SL1344/coloc/coloc_lists/AZ.coloc.txt", sep = "\t", quote = FALSE, row.names = FALSE)






#Make plots of the stronges hits
coloc_hits = dplyr::filter(coloc_res, PP.H4.abf > 0.5)

data_res = purrr::by_row(coloc_hits, ~colocMolecularQTLs(.,qtl_summary_path = qtlResults()$rna_fastqtl$naive, 
              gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz")
              ,GRCh37_variants, GRCh38_variants, qtl_type = "fastqtl")$data)
data_res = dplyr::left_join(data_res, gene_name_map)
plot_list = purrr::map(data_res$.out, ~makeColocPlot(.))
names(plot_list) = data_res$gene_name
savePlotList(plot_list, "results/SL1344/coloc/UC")



#Some example code to debug coloc scripts
gwas_top_hits = importGWASSummary(paste0(gwas_prefix, ".sorted.txt.gz"))

test_res = colocMolecularQTLsByRow(qtl_df_list$naive[1:10,], qtl_summary_path = qtlResults()$rna_fastqtl$naive, 
              gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants = GRCh37_variants, 
              GRCh38_variants = GRCh38_variants, qtl_type = "fastqtl", N_qtl = 84, cis_dist = 1e5)

data_res = purrr::by_row(qtl_df_list$naive[1:10,], ~colocMolecularQTLs(.,qtl_summary_path = qtlResults()$rna_fastqtl$naive, 
                gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz")
                ,GRCh37_variants, GRCh38_variants, qtl_type = "fastqtl")$data)
  

