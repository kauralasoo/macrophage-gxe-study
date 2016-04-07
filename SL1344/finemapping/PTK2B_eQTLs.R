library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
library("MatrixEQTL")
load_all("macrophage-gxe-study/housekeeping/")
library("ggplot2")
load_all("~/software/rasqual/rasqualTools/")
library("purrr")

#### Import data ####
#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
#Load p-values from disk
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
rasqual_min_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_fdr < 0.1)})
min_pvalue_df = plyr::ldply(rasqual_min_hits, .id = "condition_name") %>% dplyr::arrange(p_nominal)
joint_pairs = dplyr::select(min_pvalue_df, gene_id, snp_id) %>% unique() 

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Calculate R2
genotypes = vcf_file$genotypes[unique(joint_pairs$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Prune SNPs
filtered_pairs = filterHitsR2(joint_pairs, filtered_vcf$genotypes, .8)

#Find unique eQTLs for this gene:
independent_qtls = dplyr::filter(filtered_pairs, gene_id == "ENSG00000120899")
independent_qtl_res = dplyr::semi_join(min_pvalue_df, independent_qtls)

#Make effect size plots for these two qtls
naive_eQTL = plotEQTL("ENSG00000120899", "rs6987305", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
ggsave("results/SL1344/eQTLs/example_loci/PTK2B/PTK2B_naive_eQTL.pdf", naive_eQTL, width = 7, height = 7)

ifng_sl1344_eQTL = plotEQTL("ENSG00000120899", "rs1429938", combined_expression_data$cqn, vcf_file$genotypes, 
                      combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
ggsave("results/SL1344/eQTLs/example_loci/PTK2B/PTK2B_IFNg_SL1344_eQTL.pdf", ifng_sl1344_eQTL, width = 7, height = 7)


#Import pvalues for PTK2B gene
naive_pvalues = tabixFetchGenesQuick("ENSG00000120899", "databases/SL1344/naive_500kb.sorted.txt.gz", combined_expression_data$gene_metadata, cis_window = 0.5e5)[[1]] %>% 
  dplyr::mutate(condition = "naive") %>% addAssociationPosterior(n = 69)
IFNg_SL1344_pvalues = tabixFetchGenesQuick("ENSG00000120899", "databases/SL1344/IFNg_SL1344_500kb.sorted.txt.gz", combined_expression_data$gene_metadata, cis_window = 0.5e5)[[1]]  %>% 
  dplyr::mutate(condition = "IFNg_SL1344") %>% addAssociationPosterior(n = 69)
joint_pvalues = rbind(naive_pvalues, IFNg_SL1344_pvalues)

#Import Alzheimer's GWAS summary stats
igap_gwas_results = readRDS("annotations/IGAP_stage_1_2_combined_stats.rds")
gwas_pvalues = dplyr::filter(igap_gwas_results, chr == 8, pos > min(joint_pvalues$pos), pos < max(joint_pvalues$pos)) %>%
  dplyr::transmute(pos, p_nominal = igap_pvalue, condition = "AD GWAS")

#Plot all pvalues
all_pvalues = rbind(dplyr::select(joint_pvalues, pos, p_nominal, condition), gwas_pvalues)
all_pvalues$condition = factor(all_pvalues$condition, levels = c("AD GWAS", "naive", "IFNg_SL1344"))
eQTL_manhattan_plot = ggplot(all_pvalues, aes(x = pos, y = -log(p_nominal, 10))) + geom_point() + facet_grid(condition~., scales = "free_y")
ggsave("results/SL1344/eQTLs/example_loci/PTK2B/eQTL_vs_GWAS_manhattan_plot.pdf", plot = eQTL_manhattan_plot, width = 7, height = 8)

#Find the most associated ATAC peaks for the lead variants
#Construct GRanges object of SNP positions
atac_tabix_list = list(naive = "../macrophage-chromatin/results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz",
                       IFNg = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz",
                       SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz",
                       IFNg_SL1344 = "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz")
peak_midpoints = dplyr::transmute(atac_list$gene_metadata, gene_id, midpoint = end -((end-start)/2)) %>% tbl_df()

#Fetch ATAC SNPs
atac_snp_results = purrr::map(atac_tabix_list, ~tabixFetchSNPsQuick(c("rs1429938","rs6987305"),.,vcf_file$snpspos)) %>%
  ldply(.id = "condition") %>%
  dplyr::left_join(peak_midpoints)

naive_filtered = dplyr::filter(atac_snp_results, snp_id == "rs6987305")
assoc_atac_peaks = ggplot(naive_filtered, aes(x = midpoint, y = -log(p_nominal, 10))) + 
  geom_point() + 
  facet_grid(condition~.) + 
  geom_vline(x = 27350609) +
  scale_x_continuous(limits = c(min(joint_pvalues$pos),max(joint_pvalues$pos)))
ggsave("results/SL1344/eQTLs/example_loci/PTK2B/naive_lead_ATAC_peaks.pdf", plot = assoc_atac_peaks, width = 7, height = 7)

stimulated_filtered = dplyr::filter(atac_snp_results, snp_id == "rs1429938")
assoc_atac_peaks = ggplot(stimulated_filtered, aes(x = midpoint, y = -log(p_nominal, 10))) + 
  geom_point() + 
  facet_grid(condition~.) + 
  geom_vline(x = 27379689) +
  scale_x_continuous(limits = c(min(joint_pvalues$pos),max(joint_pvalues$pos)))
ggsave("results/SL1344/eQTLs/example_loci/PTK2B/stimulated_lead_ATAC_peaks.pdf", plot = assoc_atac_peaks, width = 7, height = 7)


#Make effect size plots for the ATAC QTLs
#Make effect size plots for these two qtls
naive_QTL = plotEQTL("ATAC_peak_261927", "rs6987305", atac_list$cqn, vcf_file$genotypes, 
                      atac_list$sample_metadata, atac_list$gene_metadata)
ggsave("results/SL1344/eQTLs/example_loci/PTK2B/PTK2B_enhancer1_QTL.pdf", naive_QTL, width = 7, height = 7)

ifng_sl1344_QTL = plotEQTL("ATAC_peak_261942", "rs1429938", atac_list$cqn, vcf_file$genotypes, 
                           atac_list$sample_metadata, atac_list$gene_metadata)
ggsave("results/SL1344/eQTLs/example_loci/PTK2B/PTK2B_enhancer2_eQTL.pdf", ifng_sl1344_QTL, width = 7, height = 7)


#Fetch rasqual results for the ATAC peaks
naive_atac_pvalues = tabixFetchGenesQuick("ATAC_peak_261927", "../macrophage-chromatin/results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz", atac_list$gene_metadata, cis_window = 0.5e5)[[1]] %>% 
  dplyr::mutate(condition = "naive")
IFNg_SL1344_atac_pvalues = tabixFetchGenesQuick("ATAC_peak_261942", "../macrophage-chromatin/results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz", atac_list$gene_metadata, cis_window = 0.5e5)[[1]]  %>% 
  dplyr::mutate(condition = "IFNg_SL1344")


#Test for colocalisation
a = testColoc(naive_pvalues, naive_atac_pvalues, n1 = 69, n2 = 42)
b = testColoc(naive_pvalues, IFNg_SL1344_atac_pvalues, n1 = 69, n2 = 42)
c = testColoc(IFNg_SL1344_pvalues, IFNg_SL1344_atac_pvalues, n1 = 69, n2 = 31)
c = testColoc(IFNg_SL1344_pvalues, naive_pvalues, n1 = 69, n2 = 69)



#Look at some other variants
EPHA1_pvalues = tabixFetchGenesQuick("ENSG00000229153", "databases/SL1344/naive_500kb.sorted.txt.gz", combined_expression_data$gene_metadata, cis_window = 5e5)[[1]] %>%
  dplyr::transmute(CHR = chr, BP = pos, SNP = snp_id, P = p_nominal)
write.table(EPHA1_pvalues, "EPHA1-AS.gwas", sep ="\t", quote = FALSE, row.names = FALSE)

pvalues = tabixFetchGenesQuick("ENSG00000029725", "databases/SL1344/naive_500kb.sorted.txt.gz", combined_expression_data$gene_metadata, cis_window = 5e5)[[1]] %>%
  dplyr::transmute(CHR = chr, BP = pos, SNP = snp_id, P = p_nominal)
write.table(pvalues, "RABEP.gwas", sep ="\t", quote = FALSE, row.names = FALSE)

pvalues = tabixFetchGenesQuick("ENSG00000030066", "databases/SL1344/naive_500kb.sorted.txt.gz", combined_expression_data$gene_metadata, cis_window = 5e5)[[1]] %>%
  dplyr::transmute(CHR = chr, BP = pos, SNP = snp_id, P = p_nominal)
write.table(pvalues, "NUP160.gwas", sep ="\t", quote = FALSE, row.names = FALSE)

pvalues = tabixFetchGenesQuick("ENSG00000105383", "databases/SL1344/naive_500kb.sorted.txt.gz", combined_expression_data$gene_metadata, cis_window = 5e5)[[1]] %>%
  dplyr::transmute(CHR = chr, BP = pos, SNP = snp_id, P = p_nominal)
write.table(pvalues, "CD33.gwas", sep ="\t", quote = FALSE, row.names = FALSE)


#Perform finemapping with ATAC data

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")
min_pvalues_list = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalues_hits = lapply(min_pvalues_list, function(x){dplyr::filter(x, p_fdr < 0.1)})

#Overlap
naive_credible_set = constructCredibleSet(naive_pvalues, threshold = 0.99)
naive_cs_annotated = annotateCredibleSet(naive_credible_set, atac_list$gene_metadata, atac_tabix_list$naive)

ifng_sl1344_cs = constructCredibleSet(IFNg_SL1344_pvalues, threshold = 0.99)
ifng_sl1344_cs_annotated = annotateCredibleSet(ifng_sl1344_cs, atac_list$gene_metadata, atac_tabix_list$IFNg_SL1344)




