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
#opt = list(g = "IBD", w = "2e5", p = "full", d = "../../datasets/Inflammatory_GWAS/", o = "results/Fairfax/coloc/coloc_lists/")

#Extract parameters for CMD options
gwas_id = opt$g
cis_window = as.numeric(opt$w)
phenotype = opt$p
gwas_dir = opt$d
outdir = opt$o

#Import variant information
GRCh38_variants = importVariantInformation("processed/Fairfax/merged_genotypes/fairfax_genotypes.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("processed/Fairfax/merged_genotypes/fairfax_genotypes.variant_information.txt.gz")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name","type"))

#Make a list of sample sizes
sample_sizes_all = list(full = list(CD14 = 414, IFN = 367, LPS2 = 261, LPS24 = 322),
                        shared = list(CD14 = 228, IFN = 228, LPS2 = 228, LPS24 = 228),
                        shared_84 = list(CD14 = 84, IFN = 84, LPS2 = 84, LPS24 = 84))

#Construct qtl_list
phenotype_list = list(
  full = constructQtlListForColoc(phenotype = "full", qtl_root = "processed/Fairfax/qtltools/output/", 
                                  sample_size_list = sample_sizes_all[["full"]]),
  shared = constructQtlListForColoc(phenotype = "shared", qtl_root = "processed/Fairfax/qtltools/output/", 
                                  sample_size_list = sample_sizes_all[["shared"]]),
  shared_84 = constructQtlListForColoc(phenotype = "shared_84", qtl_root = "processed/Fairfax/qtltools/output/", 
                                  sample_size_list = sample_sizes_all[["shared_84"]])
)

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
                                                      N_qtl = .y, cis_dist = cis_window))

#Export results
coloc_hits = purrr::map_df(coloc_res_list, identity, .id = "condition_name") %>% dplyr::arrange(gwas_lead)
coloc_output = file.path(outdir, paste(gwas_id, phenotype, opt$w, "txt", sep = "."))
write.table(coloc_hits, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)

#Debugging example
#colocMolecularQTLs(qtl_pairs[1,], qtl_summary_list$Ctrl, gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants, GRCh38_variants, N_qtl = 2e5)
