library("devtools")
library("plyr")
library("dplyr")
library("purrr")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
load_all("~/software/rasqual/rasqualTools/")

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

#Import eQTL lists
rna_fastqtl_100kb = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues_100kb.rds")
naive_eqtls = dplyr::filter(rna_fastqtl_100kb$naive, p_nominal < 1e-8) %>% dplyr::select(gene_id)
naive_eqtl_ranges = constructGeneRanges(naive_eqtls, tss_coords, 1e5)
rna_pvals = fastqtlTabixFetchGenes(naive_eqtl_ranges, qtlResults()$rna_fastqtl$naive)

#Convert p-values to z-scores and compute LD matrices
rna_z = purrr::map(rna_pvals, pToZ)
rna_z = rna_z[unlist(map(rna_z, nrow)) > 0] #Remove genes with 0 variants
rna_ld = purrr::map(rna_z, zToLD, vcf_file$genotypes)
saveFinemapMatrices(rna_z, "results/SL1344/finemap/", file_suffix = "z")
saveFinemapMatrices(rna_ld, "results/SL1344/finemap/", file_suffix = "ld")

#Make a master data frame pointing to all of the files
rna_meta = finemapConstructMeta(names(rna_z), "results/SL1344/finemap/", 84)
write.table(rna_meta, "results/SL1344/finemap/data.txt", sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Import eQTL lists (42 samples)
rna_fastqtl_100kb = readRDS("results/SL1344/eQTLs/fastqtl_min_pvalues_100kb_42.rds")
naive_eqtls = dplyr::filter(rna_fastqtl_100kb$naive, p_nominal < 1e-8) %>% dplyr::select(gene_id)
naive_eqtl_ranges = constructGeneRanges(naive_eqtls, tss_coords, 1e5)
rna_pvals = fastqtlTabixFetchGenes(naive_eqtl_ranges, "/Volumes/JetDrive/databases/SL1344/fastqtl/tss_42/naive_500kb_pvalues.sorted.txt.gz")

#Convert p-values to z-scores and compute LD matrices
rna_z = purrr::map(rna_pvals, pToZ)
rna_z = rna_z[unlist(map(rna_z, nrow)) > 0] #Remove genes with 0 variants
rna_ld = purrr::map(rna_z, zToLD, vcf_file$genotypes)
saveFinemapMatrices(rna_z, "results/SL1344/finemap_n42/", file_suffix = "z")
saveFinemapMatrices(rna_ld, "results/SL1344/finemap_n42/", file_suffix = "ld")

#Make a master data frame pointing to all of the files
rna_meta = finemapConstructMeta(names(rna_z), "results/SL1344/finemap_n42/", 84)
write.table(rna_meta, "results/SL1344/finemap_n42/data.txt", sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)


#Process caQTLs
atac_fastqtl_100kb = readRDS("results/ATAC/QTLs/fastqtl_min_pvalues_100kb.rds")
naive_caqtls = dplyr::filter(atac_fastqtl_100kb$naive, p_nominal < 1e-8) %>% dplyr::select(gene_id)
naive_caqtls = naive_caqtls[sample(3605,849),]
naive_caqtl_ranges = constructGeneRanges(naive_caqtls, peak_centres, 1e5)
atac_pvals = fastqtlTabixFetchGenes(naive_caqtl_ranges, qtlResults()$atac_fastqtl$naive)

#Convert p-values to z-scores and compute LD matrices
atac_z = purrr::map(atac_pvals, pToZ)
atac_z = atac_z[unlist(map(atac_z, nrow)) > 0] #Remove genes with 0 variants
atac_ld = purrr::map(atac_z, zToLD, vcf_file$genotypes)
saveFinemapMatrices(atac_z, "results/ATAC/finemap/", file_suffix = "z")
saveFinemapMatrices(atac_ld, "results/ATAC/finemap/", file_suffix = "ld")

#Construct metadata
atac_meta = finemapConstructMeta(names(atac_z), "results/ATAC/finemap/", 42)
write.table(atac_meta, "results/ATAC/finemap/data.txt", sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)


#Import FINEMAP posterior probabilities
atac_post = read.table("ATAC_posteriors.txt", stringsAsFactors = FALSE) %>% dplyr::mutate(type = "caQTLs")
rna_post = read.table("RNA_posteriors.txt", stringsAsFactors = FALSE) %>% dplyr::mutate(type = "eQTLs")
post_df = dplyr::bind_rows(atac_post, rna_post) %>% dplyr::transmute(p_more = 1-V4, type)
posterior_plot = ggplot(post_df, aes(x = p_more, color = type)) + 
  geom_density() + 
  theme_light() +
  xlab("Posterior probability")
ggsave("figures/main_figures/causal_variant_count_posterior.pdf", plot = posterior_plot, width = 5, height = 4)



