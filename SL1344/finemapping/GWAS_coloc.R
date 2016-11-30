library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("readr")
library("devtools")
library("optparse")
load_all("../seqUtils/")

#Parse command-line options
option_list <- list(
  make_option(c("-t", "--type"), type="character", default=NULL,
              help="Type of the QTLs used for coloc.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Parse variables
gwas_id = NULL
trait = opt$t
print(trait)

####### Get GWAS id form STDIN ######
f <- file("stdin")
open(f)
gwas_id = readLines(f)
close(f)
####### END #######

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name"))

#Set intput and output files
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path("databases/GWAS/summary", gwas_file_name)
coloc_output = file.path("results/SL1344/coloc/coloc_lists/", paste0(gwas_id, ".", trait, ".coloc.txt"))

#Import old and new variant coordinates
GRCh38_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("genotypes/SL1344/imputed_20151005/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz")

#Import eQTLs
qtl_min_pvalues = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues.rds")

#Prefilter potential GWAS overlaps
qtl_df_list = prefilterColocCandidates(qtl_min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.1, overlap_dist = 1e5, gwas_thresh = 1e-5)
#Identify all unique pairs
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()

#Specify location of the QTL summary files
if(trait == "eQTL"){
  qtl_summary_list = list(naive = "results/SL1344/fastqtl/output_start_end/naive_500kb_pvalues.sorted.txt.gz",
                          IFNg = "results/SL1344/fastqtl/output_start_end/IFNg_500kb_pvalues.sorted.txt.gz",
                          SL1344 = "results/SL1344/fastqtl/output_start_end/SL1344_500kb_pvalues.sorted.txt.gz",
                          IFNg_SL1344 = "results/SL1344/fastqtl/output_start_end/IFNg_SL1344_500kb_pvalues.sorted.txt.gz")
  sample_sizes = list(naive = 84, IFNg = 84, SL1344 = 84, IFNg_SL1344 = 84)
} else if (trait == "caQTL"){
  qtl_summary_list = list(naive = "results/ATAC/fastqtl/output_start_end/naive_100kb_pvalues.sorted.txt.gz",
                          IFNg = "results/ATAC/fastqtl/output_start_end/IFNg_100kb_pvalues.sorted.txt.gz",
                          SL1344 = "results/ATAC/fastqtl/output_start_end/SL1344_100kb_pvalues.sorted.txt.gz",
                          IFNg_SL1344 = "results/ATAC/fastqtl/output_start_end/IFNg_SL1344_100kb_pvalues.sorted.txt.gz")
  sample_sizes = list(naive = 42, IFNg = 41, SL1344 = 31, IFNg_SL1344 = 31)
}


#Test for coloc
coloc_res_list = purrr::map2(qtl_summary_list, sample_sizes, ~colocMolecularQTLsByRow(qtl_pairs, qtl_summary_path = .x, 
                   gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants = GRCh37_variants, 
                   GRCh38_variants = GRCh38_variants, qtl_type = "fastqtl", N_qtl = .y, cis_dist = 1e5))

coloc_hits = purrr::map_df(coloc_res_list, identity, .id = "condition_name") %>% arrange(gwas_lead)
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
  

