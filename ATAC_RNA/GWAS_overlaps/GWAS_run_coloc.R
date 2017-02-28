library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("readr")
library("devtools")
library("optparse")
load_all("../../seqUtils/")

#Parse command-line options
option_list <- list(
  make_option(c("-p", "--phenotype"), type="character", default=NULL,
              help="Type of QTLs used for coloc.", metavar = "type"),
  make_option(c("-w", "--window"), type="character", default=NULL,
              help="Size of the cis window.", metavar = "type"),
  make_option(c("-g", "--gwas"), type="character", default=NULL,
              help="Name of the GWAS trait", metavar = "type"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Path to GWAS summary stats directory.", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Debugging
#opt = list(g = "IBD", w = "2e5", p = "SL1344", d = "databases/GWAS/summary", o = "results/acLDL/coloc/coloc_lists/")

#Extract parameters for CMD options
gwas_id = opt$g
cis_window = as.numeric(opt$w)
phenotype = opt$p
gwas_dir = opt$d
outdir = opt$o

#Specify list of phenotypes
phenotype_list = list(
  SL1344 = list(
    min_pvalues = readRDS("../results/SL1344/eQTLs/fastqtl_min_pvalues.rds") %>%
      purrr::map(~dplyr::transmute(., phenotype_id = gene_id, snp_id, p_fdr)),
    qtl_summary_list = list(naive = "../results/SL1344/fastqtl/output_start_end/naive_500kb_pvalues.sorted.txt.gz",
                            IFNg = "../results/SL1344/fastqtl/output_start_end/IFNg_500kb_pvalues.sorted.txt.gz",
                            SL1344 = "../results/SL1344/fastqtl/output_start_end/SL1344_500kb_pvalues.sorted.txt.gz",
                            IFNg_SL1344 = "../results/SL1344/fastqtl/output_start_end/IFNg_SL1344_500kb_pvalues.sorted.txt.gz"),
    sample_sizes = list(naive = 84, IFNg = 84, SL1344 = 84, IFNg_SL1344 = 84),
    QTLTools = FALSE
  ),
  ATAC = list(
    min_pvalues = readRDS("../results/ATAC/QTLs/fastqtl_min_pvalues.rds"),
    qtl_summary_list = list(naive = "../results/ATAC/fastqtl/output/naive_500kb_pvalues.sorted.txt.gz",
                          IFNg = "../results/ATAC/fastqtl/output/IFNg_500kb_pvalues.sorted.txt.gz",
                          SL1344 = "../results/ATAC/fastqtl/output/SL1344_500kb_pvalues.sorted.txt.gz",
                          IFNg_SL1344 = "../results/ATAC/fastqtl/output/IFNg_SL1344_500kb_pvalues.sorted.txt.gz"),
  sample_sizes = list(naive = 42, IFNg = 41, SL1344 = 31, IFNg_SL1344 = 31),
  QTLTools = FALSE)
)

#Import variant information
GRCh38_variants = importVariantInformation("../genotypes/acLDL/imputed_20151005/imputed.70_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("../genotypes/acLDL/imputed_20151005/GRCh37/imputed.70_samples.variant_information.GRCh37.txt.gz")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name"))

#Extract the phenotype of interest
phenotype_values = phenotype_list[[phenotype]]

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path(gwas_dir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = prefilterColocCandidates(phenotype_values$min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.1, 
                                       overlap_dist = 1e5, gwas_thresh = 1e-5)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()

#Test for coloc
coloc_res_list = purrr::map2(phenotype_values$qtl_summary_list, phenotype_values$sample_sizes, 
                             ~colocMolecularQTLsByRow(qtl_pairs, qtl_summary_path = .x, 
                                                      gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), 
                                                      GRCh37_variants = GRCh37_variants,
                                                      GRCh38_variants = GRCh38_variants, 
                                                      N_qtl = .y, cis_dist = cis_window, 
                                                      QTLTools = phenotype_values$QTLTools))

#Export results
coloc_hits = purrr::map_df(coloc_res_list, identity, .id = "condition_name") %>% dplyr::arrange(gwas_lead)
coloc_output = file.path(outdir, paste(gwas_id, phenotype, opt$w, "txt", sep = "."))
write.table(coloc_hits, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)

#Make plots of the stronges hits
#coloc_hits = dplyr::filter(coloc_res, PP.H4.abf > 0.5)

#data_res = purrr::by_row(coloc_hits, ~colocMolecularQTLs(.,qtl_summary_path = qtlResults()$rna_fastqtl$naive, 
#              gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz")
#              ,GRCh37_variants, GRCh38_variants, qtl_type = "fastqtl")$data)
#data_res = dplyr::left_join(data_res, gene_name_map)
#plot_list = purrr::map(data_res$.out, ~makeColocPlot(.))
#names(plot_list) = data_res$gene_name
#savePlotList(plot_list, "results/SL1344/coloc/UC")



#Some example code to debug coloc scripts
#gwas_top_hits = importGWASSummary(paste0(gwas_prefix, ".sorted.txt.gz"))

#test_res = colocMolecularQTLsByRow(qtl_df_list$naive[1:10,], qtl_summary_path = qtlResults()$rna_fastqtl$naive, 
#              gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants = GRCh37_variants, 
#              GRCh38_variants = GRCh38_variants, qtl_type = "fastqtl", N_qtl = 84, cis_dist = 1e5)

#data_res = purrr::by_row(qtl_df_list$naive[1:10,], ~colocMolecularQTLs(.,qtl_summary_path = qtlResults()$rna_fastqtl$naive, 
#                gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz")
#                ,GRCh37_variants, GRCh38_variants, qtl_type = "fastqtl")$data)
  

