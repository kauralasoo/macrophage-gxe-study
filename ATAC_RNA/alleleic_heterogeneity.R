library("devtools")
library("plyr")
library("dplyr")
library("purrr")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")

#Helper functions
pToZ <- function(pvalue_df, n_top_snps = 200){
  result = pvalue_df %>% 
    dplyr::arrange(p_nominal) %>%
    head(n_top_snps) %>% 
    dplyr::arrange(pos) %>%
    dplyr::mutate(z = qnorm(p_nominal)*sign(beta)) %>%
    dplyr::select(snp_id, z)
  return(result)
}

zToLD <- function(z_df, genotypes){
  cor_mat = cor(t(genotypes[z_df$snp_id,]), method = "pearson", use = "pairwise.complete.obs")
  return(cor_mat)
}

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
atac_list = readRDS("results/ATAC/ATAC_combined_accessibility_data.rds")

#Calculate coordinates of the TSSs
tss_coords = dplyr::mutate(combined_expression_data$gene_metadata, tss = ifelse(strand == 1, start, end)) %>%
  dplyr::mutate(start = tss, end = tss)
peak_centres = dplyr::mutate(atac_list$gene_metadata, centre = floor((start + end)/2)) %>% 
  dplyr::mutate(start = centre, end = centre)

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import variant information
snp_info = importVariantInformation("genotypes/SL1344/imputed_20151005/imputed.86_samples.variant_information.txt.gz")

#Import QTL lists
rna_fastqtl_100kb = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues_100kb.rds")
naive_eqtls = dplyr::filter(rna_fastqtl_100kb$naive, p_nominal < 1e-8) %>% dplyr::select(gene_id)
naive_eqtl_ranges = constructGeneRanges(naive_eqtls, tss_coords, 1e5)
rna_pvals = fastqtlTabixFetchGenes(naive_eqtl_ranges, qtlResults()$rna_fastqtl$naive)

#Convert p-values to z-scores and compute LD matrices
rna_z = purrr::map(rna_pvals, pToZ)
rna_z = rna_z[unlist(map(rna_z, nrow)) > 0] #Remove genes with 0 variants
rna_ld = purrr::map(rna_z, zToLD, vcf_file$genotypes)



#Construct gene ranges
ptk2b_region = constructGeneRanges(data_frame(gene_id = "ENSG00000120899"), tss_coords, 2e5)
peak1_region = constructGeneRanges(data_frame(gene_id = "ATAC_peak_261927"), peak_centres, 1e5)
ptk2b_region = constructGeneRanges(data_frame(gene_id = "ENSG00000170458"), tss_coords, 2e5)

rna_variants = fastqtlTabixFetchGenes(ptk2b_region, qtlResults()$rna_fastqtl$naive)[[1]] %>% 
  dplyr::arrange(p_nominal) %>%
  head(200) %>% 
  dplyr::arrange(pos) %>%
  dplyr::mutate(z = qnorm(p_nominal)*sign(beta))
write.table(dplyr::select(rna_variants, snp_id, z), "../../software/finemap_v1.1_MacOSX/example/rna.z", sep = " ", 
            quote = FALSE, row.names = FALSE, col.names = FALSE) 

#Calculate correlations between variants
rna_cor = cor(t(vcf_file$genotypes[rna_variants$snp_id,]), method = "pearson")
write.table(rna_cor, "../../software/finemap_v1.1_MacOSX/example/rna.ld", 
            quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

atac_variants = fastqtlTabixFetchGenes(peak1_region, qtlResults()$atac_fastqtl$naive)[[1]] %>% 
  dplyr::arrange(p_nominal) %>%
  head(100) %>%
  dplyr::arrange(pos) %>%
  dplyr::mutate(z = qnorm(p_nominal)*sign(beta))
write.table(dplyr::select(atac_variants, snp_id, z), "../../software/finemap_v1.1_MacOSX/example/atac.z", sep = " ", 
            quote = FALSE, row.names = FALSE, col.names = FALSE) 

atac_cor = cor(t(vcf_file$genotypes[atac_variants$snp_id,]), method = "pearson")
write.table(atac_cor, "../../software/finemap_v1.1_MacOSX/example/atac.ld", 
            quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

