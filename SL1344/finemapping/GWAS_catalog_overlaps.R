library("purrr")
library("devtools")
library("plyr")
library("dplyr")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")

#Load the raw eQTL dataset
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
                                                                 levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#Import QTL variants
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
variable_qtls = readRDS("results/SL1344/eQTLs/appeat_disappear_eQTLs.rds")
rasqual_qtl_df = extractQTLsFromList(rasqual_min_pvalues, fdr_cutoff = 0.1)

#Extract pairs of SNPs
joint_pairs = dplyr::select(rasqual_qtl_df, gene_id, snp_id) %>% unique() 
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

#### GWAS overlaps ####
#Import GWAS catalog
filtered_catalog = readRDS("annotations/gwas_catalog_v1.0.1-downloaded_2016-03-02.filtered.rds")
n_loci = dplyr::group_by(filtered_catalog, trait) %>% dplyr::summarise(n_loci = length(catalog_snp_id))

#All GWAS overlaps
all_olaps = findGWASOverlaps(filtered_pairs, filtered_catalog, vcf_file, min_r2 = 0.7)
all_gwas_hits = dplyr::left_join(all_olaps, gene_name_map, by = "gene_id") %>%
  dplyr::select(gene_name, gene_id, snp_id, gwas_snp_id, R2, trait, gwas_pvalue)
write.table(all_gwas_hits, "results/SL1344/eQTLs/all_gwas_overlaps.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Rank traits by overlap size
ranked_traits = rankTraitsByOverlapSize(dplyr::filter(all_gwas_hits, R2 > 0.78), filtered_catalog, min_overlap = 5)
write.table(ranked_traits, "results/SL1344/eQTLs/relative_gwas_overlaps_R08.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Find which QTLs are condition specific
appear_eqtls = dplyr::select(variable_qtls$appear, gene_id, snp_id, cluster_id) %>% unique()
appear_gwas_hits = dplyr::semi_join(all_gwas_hits, appear_eqtls, by = c("gene_id"))

#Explore IBD traits
ibd_appear_olaps = dplyr::filter(appear_gwas_hits, trait %in% c("Ulcerative colitis", "Inflammatory bowel disease","Crohn's disease"))

#Import GWAS summary staistics around specific loci
snp_ranges = dplyr::filter(vcf_file$snpspos, snpid %in% unique(ibd_appear_olaps$snp_id)) %>% 
  dplyr::transmute(seqnames = chr, start = pos - 100000, end = pos + 100000, strand = "+", snp_id = snpid) %>% 
  dataFrameToGRanges()
ibd_pvalues = scanTabixDataFrame("/Volumes/JetDrive/databases/GWAS/IBD/EUR.IBD.gwas.bed.GRCh38.sorted.txt.gz", snp_ranges, col_names = FALSE)









#Are interaction QTLs more enriched for certain traits
immune_traits = c("Rheumatoid arthritis","Inflammatory bowel disease","Crohn's disease",
                  "Ulcerative colitis","Type 2 diabetes", "Inflammatory skin disease",
                  "Psoriasis","Systemic lupus erythematosus","Celiac disease","Type 1 diabetes",
                  "Primary biliary cirrhosis","Atopic dermatitis","Systemic sclerosis",
                  "Ankylosing spondylitis","Lupus nephritis in systemic lupus erythematosus",
                  "Systemic lupus erythematosus and Systemic sclerosis",
                  "Celiac disease or Rheumatoid arthritis","Amyotrophic lateral sclerosis (age of onset)",
                  "Type 1 diabetes and autoimmune thyroid diseases")
secondary_immmune_traits = c("Allergic rhinitis","Asthma","Immune response to smallpox vaccine (IL-6)
                             ", "Type 1 diabetes autoantibodies","C-reactive protein levels","Response to tocilizumab in rheumatoid arthritis
                             ","Immunoglobulin A", "C-reactive protein","Monocyte count")

selected_immune_traits = c("Rheumatoid arthritis","Inflammatory bowel disease","Crohn's disease",
                           "Ulcerative colitis","Type 2 diabetes")

brain_traits = c("Schizophrenia","Parkinson's disease","Alzheimer's disease (late onset)","Alzheimer's disease","Psychosis and Alzheimer's disease","Alzheimer's disease (age of onset)")

lipid_traits = c("Cholesterol, total", "HDL cholesterol", "LDL cholesterol","Coronary heart disease",
                 "Triglycerides","Trans fatty acid levels","Cardiovascular disease risk factors",
                 "Very long-chain saturated fatty acid levels (fatty acid 20:0)",
                 "Lipoprotein-associated phospholipase A2 activity and mass",
                 "Presence of antiphospholipid antibodies")

dplyr::filter(all_gwas_hits, trait %in% immune_traits, R2 > 0.8)$snp_id %>% unique() %>% length()
dplyr::filter(interaction_gwas_hits, trait %in% immune_traits, R2 > 0.8)$snp_id %>% unique() %>% length()
dplyr::filter(all_gwas_hits, trait %in% immune_traits, R2 > 0.8, snp_id %in% appear_qtls$snp_id)$snp_id %>% unique() %>% length()

#Immune
17/1233
57/5660
fisher.test(t(matrix(c(17,1233,40,5660-1233),nrow=2)))

dplyr::filter(all_gwas_hits, trait %in% selected_immune_traits, R2 > 0.8)$gwas_snp_id %>% unique() %>% length()
dplyr::filter(interaction_gwas_hits, trait %in% selected_immune_traits, R2 > 0.8)$gwas_snp_id %>% unique() %>% length()
dplyr::filter(all_gwas_hits, trait %in% selected_immune_traits, R2 > 0.8, snp_id %in% appear_qtls$snp_id)$snp_id %>% unique() %>% length()

#Immune subset
14/1233
37/5660
fisher.test(t(matrix(c(14,1233,23,5660-1233),nrow=2)))


dplyr::filter(all_gwas_hits, trait %in% brain_traits, R2 > 0.8)$snp_id %>% unique() %>% length()
dplyr::filter(interaction_gwas_hits, trait %in% brain_traits, R2 > 0.8)$snp_id %>% unique() %>% length()

#Brain
6/1233
16/5660
fisher.test(t(matrix(c(6,1233,10,5660-1233),nrow=2)))


dplyr::filter(all_gwas_hits, trait %in% lipid_traits, R2 > 0.8)$snp_id %>% unique() %>% length()
dplyr::filter(interaction_gwas_hits, trait %in% lipid_traits, R2 > 0.8)$snp_id %>% unique() %>% length()


