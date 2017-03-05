library("devtools")
library("plyr")
library("dplyr")
library("purrr")
library("ggplot2")
load_all("../../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")

#Import R2 pairs of QTLs
rna_atac_overlaps = readRDS("results/ATAC_RNA_overlaps/QTL_overlap_list_R2.rds")

#Import variant informaton
GRCh38_variants = importVariantInformation("../genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

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
    min_pvalues = readRDS("../results/ATAC/QTLs/fastqtl_min_pvalues.rds") %>%
      purrr::map(~dplyr::transmute(., phenotype_id = gene_id, snp_id, p_fdr)),
    qtl_summary_list = list(naive = "../results/ATAC/fastqtl/output/naive_500kb_pvalues.sorted.txt.gz",
                            IFNg = "../results/ATAC/fastqtl/output/IFNg_500kb_pvalues.sorted.txt.gz",
                            SL1344 = "../results/ATAC/fastqtl/output/SL1344_500kb_pvalues.sorted.txt.gz",
                            IFNg_SL1344 = "../results/ATAC/fastqtl/output/IFNg_SL1344_500kb_pvalues.sorted.txt.gz"),
    sample_sizes = list(naive = 42, IFNg = 41, SL1344 = 31, IFNg_SL1344 = 31),
    QTLTools = FALSE)
)

rna_atac_overlaps


#Set test data
gene_df = dplyr::filter(phenotype_list$SL1344$min_pvalues$IFNg, phenotype_id == "ENSG00000144228")
peaks_df = dplyr::filter(rna_atac_overlaps, gene_id == "ENSG00000144228") %>% dplyr::select(peak_id) %>% unique()

coloc_res = colocGeneAgainstPeaks(gene_df, peaks_df, qtlResults()$rna_fastqtl$IFNg, 
                                  qtlResults()$atac_fastqtl$IFNg, n_eqtl = 84,
                                  n_caqtl = 41, variant_information = GRCh38_variants)


  
  









