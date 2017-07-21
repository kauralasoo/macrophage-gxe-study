library("dplyr")
library("tidyr")
library("purrr")
library("UpSetR")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import SummarizedExperiment
se_fairfax = readRDS("results/Fairfax/expression_data.SummarizedExperiment.rds")
expression_mat = assays(se_fairfax)$exprs
sample_metadata = colData(se_fairfax) %>% tbl_df2()
gene_metadata = rowData(se_fairfax) %>% tbl_df2() %>% dplyr::mutate(gene_id = probe_id)

#Import genotypes
vcf_file = seqUtils::gdsToMatrix("processed/Fairfax/geno_by_chr/9.gds")
variant_information = importVariantInformation("processed/Fairfax/merged_genotypes/fairfax_genotypes.variant_information.txt.gz")

#Make a QTL plot
TRAF1_data = constructQtlPlotDataFrame("ILMN_1698218", "rs10985070", 
                                       expression_mat, 
                                       vcf_file$genotypes, 
                                       sample_metadata, 
                                       gene_metadata) %>%
  dplyr::left_join(constructGenotypeText("rs10985070", variant_information), by = "genotype_value")

TRAF1_qtl_plot = plotQtlCol(TRAF1_data, scales = "free_y")
ggsave("figures/supplementary/fairfax_TRAF1_boxplot.pdf", width = 3, height = 5)
  

