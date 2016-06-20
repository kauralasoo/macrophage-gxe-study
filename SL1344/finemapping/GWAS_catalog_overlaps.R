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

#Import eQTL variants
rasqual_min_pvalues = readRDS("results/SL1344/eQTLs/rasqual_min_pvalues.rds")
variable_qtls = readRDS("results/SL1344/eQTLs/appeat_disappear_eQTLs.rds")
rasqual_qtl_df = extractQTLsFromList(rasqual_min_pvalues, fdr_cutoff = 0.1)

#Extract pairs of SNPs
joint_pairs = dplyr::select(rasqual_qtl_df, gene_id, snp_id) %>% unique() 
filtered_pairs = filterHitsR2(joint_pairs, vcf_file$genotypes, .8)

#Import ATAC QTL variants
atac_min_pvalues = readRDS("../macrophage-chromatin/results/ATAC/QTLs/rasqual_min_pvalues.rds")
atac_qtl_df = extractQTLsFromList(atac_min_pvalues, fdr_cutoff = 0.1)
atac_joint_pairs = dplyr::select(atac_qtl_df, gene_id, snp_id) %>% unique() 
atac_filtered_pairs = filterHitsR2(atac_joint_pairs, vcf_file$genotypes, .8)

#### GWAS overlaps ####
#Import GWAS catalog
filtered_catalog = readRDS("annotations/gwas_catalog_v1.0.1-downloaded_2016-03-02.filtered.rds")
n_loci = dplyr::group_by(filtered_catalog, trait) %>% dplyr::summarise(n_loci = length(catalog_snp_id))

#All GWAS eQTL overlaps
all_olaps = findGWASOverlaps(filtered_pairs, filtered_catalog, vcf_file, min_r2 = 0.7)
saveRDS(all_olaps, "results/SL1344/eQTLs/GWAS_overlaps/RNA_gwas_overlaps.rds")
all_olaps = readRDS("results/SL1344/eQTLs/GWAS_overlaps/RNA_gwas_overlaps.rds")
all_gwas_hits = dplyr::left_join(all_olaps, gene_name_map, by = "gene_id") %>%
  dplyr::select(gene_name, gene_id, snp_id, gwas_snp_id, R2, trait, gwas_pvalue)
write.table(all_gwas_hits, "results/SL1344/eQTLs/GWAS_overlaps/all_gwas_overlaps.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#All GWAS caQTL overlaps
atac_olaps = findGWASOverlaps(atac_filtered_pairs, filtered_catalog, vcf_file, min_r2 = 0.7)
saveRDS(atac_olaps, "results/SL1344/eQTLs/ATAC_gwas_overlaps.rds")
atac_olaps = readRDS("results/SL1344/eQTLs/GWAS_overlaps/ATAC_gwas_overlaps.rds")
atac_gwas_hits = dplyr::transmute(atac_olaps, gene_name = gene_id, gene_id, chr = chr.x, pos = pos.x, snp_id, gwas_snp_id, R2, trait, gwas_pvalue)

#Rank traits by overlap size
ranked_traits = rankTraitsByOverlapSize(dplyr::filter(all_gwas_hits, R2 > 0.78), filtered_catalog, min_overlap = 5)
write.table(ranked_traits, "results/SL1344/eQTLs/GWAS_overlaps/relative_gwas_overlaps_R08.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Find which QTLs are condition specific
appear_eqtls = dplyr::select(variable_qtls$appear, gene_id, snp_id, cluster_id) %>% unique()
appear_gwas_hits = dplyr::semi_join(all_gwas_hits, appear_eqtls, by = c("gene_id"))

#Explore IBD traits
appear_olaps = dplyr::filter(appear_gwas_hits, trait %in% c("Ulcerative colitis", "Inflammatory bowel disease","Crohn's disease"))

#Plot QTL overlaps for IBD
#Extract overlapping IBD GWAS hits and consturct Granges
ibd_appear_olaps = dplyr::filter(appear_gwas_hits, trait %in% c("Inflammatory bowel disease")) %>%
  addVariantCoords(vcf_file$snpspos)
ibd_ranges = constructGWASRanges(ibd_appear_olaps, 200000)
gwas_colnames = c("chr","pos","pos2","snp_id","A1","A2", "INFO","OR","SE","p_nominal")
ibd_pvalues = scanTabixDataFrame("/Volumes/JetDrive/databases/GWAS/IBD/EUR.IBD.gwas.bed.GRCh38.sorted.txt.gz", ibd_ranges, col_names = gwas_colnames)
ibd_meta_pvalues = scanTabixDataFrame("/Volumes/JetDrive/databases/GWAS/IBD/meta-analysis/IBD_trans_ethnic_bed.GRCh38.sorted.txt.gz", ibd_ranges, col_names = gwas_colnames)

plot(ibd_meta_pvalues[[5]]$pos, -log(ibd_meta_pvalues[[5]]$p_nominal, 10))


#Fetch gwas hits
eqtl_pvalues = tabixFetchGenes(ibd_ranges, "results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.sorted.txt.gz")
plot(eqtl_pvalues[[5]]$pos, -log(eqtl_pvalues[[5]]$p_nominal, 10))

#Fetch fastQTL p-values
fastqtl_pvalues = fastqtlTabixFetchGenes(ibd_ranges, "/Volumes/JetDrive/databases/SL1344/fastqtl/IFNg_SL1344_500kb_pvalues.sorted.txt.gz")
plot(fastqtl_pvalues[[5]]$pos, -log(fastqtl_pvalues[[5]]$p_nominal, 10))

#Fetch atac hits
atac_ranges = dplyr::filter(atac_gwas_hits, trait == "Inflammatory bowel disease", gwas_snp_id == "rs12654812") %>% constructGWASRanges(200000)
atac_pvalues = tabixFetchGenes(atac_ranges, "../macrophage-chromatin/results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz")
atac_fastqtl_pvalues = fastqtlTabixFetchGenes(atac_ranges, "/Volumes/JetDrive/databases/ATAC/naive_100kb_pvalues.sorted.txt.gz")


#RGS14 example
eqtl = eqtl_pvalues[[5]] %>%  dplyr::mutate(phenotype = "RGS14 eQTL (RASQUAL)") %>% dplyr::select(chr, pos, phenotype, p_nominal)
fastqtl = fastqtl_pvalues[[5]] %>%  dplyr::mutate(phenotype = "RGS14 eQTL (linear)") %>% dplyr::select(chr, pos, phenotype, p_nominal)
ibd = ibd_pvalues[[5]] %>% dplyr::mutate(phenotype = "IBD") %>% dplyr::select(chr, pos, phenotype, p_nominal) %>%
  dplyr::semi_join(eqtl, by = c("chr","pos"))
ibd_meta = dplyr::transmute(ibd_meta_pvalues[[5]], chr, pos, phenotype = "IBD meta-analysis", p_nominal)
atac = dplyr::transmute(atac_fastqtl_pvalues[[1]], chr, pos, phenotype = "RGS14 caQTL", p_nominal)
res = rbind(ibd, eqtl, fastqtl, ibd_meta, atac)
rgs14_plot = ggplot(res, aes(x = pos, y = -log(p_nominal,10))) + 
  geom_point() + 
  facet_wrap(~phenotype, ncol = 1, scale = "free_y") + 
  xlab("Chromosome 5 position")
rgs14_plot
ggsave("results/SL1344/eQTLs/GWAS_overlaps/RGS14_manhattan.pdf", plot = rgs14_plot, width = 6, height = 8)



#Ulcerative colitis
uc_ranges = dplyr::filter(appear_gwas_hits, trait %in% c("Ulcerative colitis")) %>%
  addVariantCoords(vcf_file$snpspos) %>%
  constructGWASRanges(100000)
uc_pvalues = scanTabixDataFrame("/Volumes/JetDrive/databases/GWAS/IBD/EUR.UC.gwas.bed.GRCh38.sorted.txt.gz", uc_ranges, col_names = gwas_colnames)

#Chron's disease
cd_ranges = dplyr::filter(appear_gwas_hits, trait %in% c("Crohn's disease")) %>%
  addVariantCoords(vcf_file$snpspos) %>%
  constructGWASRanges(100000)
cd_pvalues = scanTabixDataFrame("/Volumes/JetDrive/databases/GWAS/IBD/EUR.CD.gwas.bed.GRCh38.sorted.txt.gz", cd_ranges, col_names = gwas_colnames)
plot(cd_pvalues[[2]]$pos, -log(cd_pvalues[[2]]$p_nominal, 10))



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


