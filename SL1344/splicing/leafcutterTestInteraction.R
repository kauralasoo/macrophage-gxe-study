library("readr")
library("dplyr")
library("tidyr")
library("limma")
library("devtools")
library("purrr")
load_all("../seqUtils/")
library("GenomicRanges")

#Import gene metadata
combined_expression_data_filtered = readRDS("results/SL1344/combined_expression_data.rds")
gene_ranges = dplyr::transmute(combined_expression_data_filtered$gene_metadata, gene_id, seqnames = chr, start, end, strand = "+") %>% 
  dataFrameToGRanges()

#Import proportion data
prop_list = readRDS("results/SL1344/combined_proportions.row_quantile.rds")
cluster_meta = dplyr::select(prop_list$gene_metadata, gene_id, cluster_id, cluster_size)

#Link leafCutter clusters to genes
cluster_ranges = tidyr::separate(cluster_meta, gene_id, c("seqnames", "start", "end", "cluster"), sep = ":") %>% 
  dplyr::select(-cluster) %>% 
  dplyr::group_by(cluster_id) %>% 
  dplyr::summarise(seqnames = seqnames[1], start = min(as.numeric(start)), end = max(as.numeric(end)), cluster_size = cluster_size[1], strand = "+") %>% 
  dataFrameToGRanges()
olaps = findOverlaps(cluster_ranges, gene_ranges)

cluster_gene_map = data_frame(cluster_id = cluster_ranges[queryHits(olaps),]$cluster_id, ensembl_gene_id = gene_ranges[subjectHits(olaps),]$gene_id) %>% 
  dplyr::group_by(cluster_id) %>% 
  dplyr::mutate(overlap_count = length(ensembl_gene_id)) %>% 
  dplyr::filter(overlap_count == 1) %>% 
  dplyr::select(-overlap_count)

cluster_named_meta = dplyr::left_join(cluster_meta, cluster_gene_map, by = "cluster_id")

#Import min p-values
fastqtl_pvalue_list = readRDS("results/SL1344/leafcutter/leafcutter_min_pvalues.rds")

#Add metadata to leafcutter QTLs
fastqtl_pvalue_meta = purrr::map(fastqtl_pvalue_list, ~dplyr::left_join(., cluster_named_meta, by = "gene_id"))

#Apply bonferroni correction for p-values within cluster
fastqtl_bonferroni = purrr::map(fastqtl_pvalue_meta, ~dplyr::group_by(., cluster_id) %>% 
             dplyr::mutate(n_transcripts = length(gene_id)) %>% 
             dplyr::arrange(cluster_id, p_beta) %>% dplyr::filter(row_number() == 1) %>% 
             dplyr::ungroup() %>% dplyr::mutate(p_bonferroni = p_beta * n_transcripts) %>% 
             dplyr::mutate(p_bonferroni = pmin(p_bonferroni, 1)) %>% 
             dplyr::mutate(p_fdr = p.adjust(p_bonferroni, method = "fdr")))
saveRDS(fastqtl_bonferroni,"results/SL1344/leafcutter/leafcutter_cluster_min_pvalues.rds")

#%>% 
#             dplyr::filter(p_fdr < 0.1))
fastqtl_bonferroni_df = purrr::map_df(fastqtl_bonferroni, identity, .id = "condition_name")

#Extract pairs
joint_pairs = dplyr::select(fastqtl_bonferroni_df, gene_id, snp_id) %>% unique()

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Filter VCF
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

#Set up formulas for interaction testing
covariate_names = c("norm_PC1", "norm_PC2", "norm_PC3","norm_PC4", "norm_PC5", "sex_binary")
formula_qtl = as.formula(paste("expression ~ genotype + condition_name ", 
                               paste(covariate_names, collapse = " + "), sep = "+ "))
formula_interaction = as.formula(paste("expression ~ genotype + condition_name + condition_name:genotype ", 
                                       paste(covariate_names, collapse = " + "), sep = "+ "))

#Test for interactions
interaction_results = testMultipleInteractions(tbl_df(filtered_pairs), prop_list$cqn, prop_list$sample_metadata, 
                                               filtered_vcf, formula_qtl, formula_interaction)
interaction_df = postProcessInteractionPvalues(interaction_results) %>% 
  dplyr::left_join(dplyr::select(cluster_meta, gene_id, cluster_id)) %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::mutate(cluster_size = length(cluster_id)) %>% 
  dplyr::mutate(p_bonferroni = p_nominal * cluster_size) %>%
  dplyr::mutate(p_bonferroni = pmin(p_bonferroni, 1)) %>%
  dplyr::mutate(p_fdr = p.adjust(p_bonferroni, method = "fdr")) %>%
  dplyr::arrange(cluster_id, p_fdr) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(p_fdr)

hist(interaction_df$p_bonferroni, breaks = 40)

dplyr::filter(interaction_df, p_fdr < 0.1)
700/1893

#Make a Q-Q plot for the interaction p-values
qq_df = dplyr::mutate(interaction_df, p_eigen = p_bonferroni) %>% 
  dplyr::arrange(p_eigen) %>% 
  addExpectedPvalue()
qq_plot = ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_bonferroni,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")
ggsave("figures/supplementary/leafcutter_interaction_qqplot.pdf", plot = qq_plot, width = 4.5, height = 4.5)

#Make a couple of plots
gene_metadata = dplyr::mutate(prop_list$gene_metadata, gene_name = gene_id)
plotEQTL("6:33085889:33086219:clu_16781", "rs34544512", prop_list$cqn, vcf_file$genotypes, 
         prop_list$sample_metadata, gene_metadata)




#Example of different variance between conditons after quantile normalisation
plotEQTL("1:179884769:179889313:clu_5474", "rs2245425", prop_list$tpm, vcf_file$genotypes, 
         prop_list$sample_metadata, gene_metadata)
plotEQTL("1:179884769:179889313:clu_5474", "rs2245425", prop_list$cqn, vcf_file$genotypes, 
         prop_list$sample_metadata, gene_metadata)

