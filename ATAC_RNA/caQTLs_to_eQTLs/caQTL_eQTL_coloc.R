library("devtools")
library("dplyr")
library("purrr")
library("coloc")
library("readr")
library("optparse")
load_all("../../seqUtils/")

#Parse command-line options
option_list <- list(
  make_option(c("-c", "--chr"), type="character", default=NULL,
              help="Chromosome used for coloc analysis.", metavar = "type"),
  make_option(c("-s", "--stimulation"), type="character", default=NULL,
              help="Stimulation condition for coloc analysis.", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))
#opt = list(c = "10", s = "IFNg_SL1344", outdir = "../processed/ATAC_RNA/eQTL_caQTL_coloc/")

#Extract options
selected_chr = opt$c
selected_condition = opt$s
outdir = opt$o

#Import R2 pairs of QTLs
rna_atac_overlaps = readRDS("../results/ATAC_RNA_overlaps/QTL_overlap_list_R2.rds")

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

#Extract parameters
rna_summaries = phenotype_list$SL1344$qtl_summary_list[[selected_condition]]
atac_summaries = phenotype_list$ATAC$qtl_summary_list[[selected_condition]]
#rna_summaries = qtlResults()$rna_fastqtl[[selected_condition]]
#atac_summaries = qtlResults()$atac_fastqtl[[selected_condition]]
n_eqtl = phenotype_list$SL1344$sample_sizes[[selected_condition]]
n_caqtl = phenotype_list$ATAC$sample_sizes[[selected_condition]]

#Run coloc over all possible pairs
res = rna_atac_overlaps %>% 
  dplyr::filter(chr == selected_chr) %>%
  #dplyr::filter(gene_id %in% c("ENSG00000120594","ENSG00000165997")) %>%
  dplyr::mutate(gene_id1 = gene_id, snp_id1 = snp_id) %>%
  group_by(gene_id1, snp_id1) %>%
  purrr::by_slice(function(slice){
    gene_df = dplyr::transmute(slice, phenotype_id = gene_id, snp_id) %>% unique()
    peak_df = dplyr::select(slice, peak_id) %>% unique()
    coloc_res = colocGeneAgainstPeaks(gene_df, peak_df, eqtl_summaries = rna_summaries, 
                                      caqtl_summaries = atac_summaries, n_eqtl = n_eqtl,
                                      n_caqtl = n_caqtl, variant_information = GRCh38_variants)
  })

#Merge results together
coloc_results = purrr::map_df(res$.out, identity)
coloc_output = file.path(outdir, paste("eQTL_caQTL_coloc", selected_condition, selected_chr, "txt", sep = "."))
write.table(coloc_results, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)

