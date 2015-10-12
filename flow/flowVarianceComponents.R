library("dplyr")
library("lme4")
library("devtools")
load_all("../seqUtils/")
library("GWASTools")
library("ggplot2")

#Import data
flow_purity = readRDS("macrophage-gxe-study/data/covariates/flow_cytometry_purity.rds")
line_medatada = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")
vcf_file = readRDS("genotypes/SL1344/array_genotypes.59_samples.imputed.vcfToMatrix.rds")
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")

#Extract gene coords
cd14_coords = dplyr::filter(eqtl_data_list$genepos, geneid == "ENSG00000170458")[,c("chr","right")]
fcgr3a_coords = dplyr::filter(eqtl_data_list$genepos, geneid == "ENSG00000203747")[,c("chr","right")]
fcgr3b_coords = dplyr::filter(eqtl_data_list$genepos, geneid == "ENSG00000162747")[,c("chr","right")]
mrc1_coords = dplyr::filter(eqtl_data_list$genepos, geneid == "ENSG00000260314")[,c("chr","left")]

#Map channel to marker
channel_marker_map = data_frame(channel = c("APC.A","PE.A","Pacific.Blue.A"), marker = c("CD206","CD16","CD14"))
unique_lines = dplyr::select(line_medatada, line_id, donor, genotype_id) %>% unique()
flow_data = dplyr::left_join(flow_purity, channel_marker_map, by = "channel") %>%
  dplyr::left_join(unique_lines, by = "donor") %>%
  dplyr::mutate(donor = ifelse(donor == "fpdj", "nibo",donor)) %>% #fpdj and nibo are the same donors
  dplyr::mutate(intensity = mean2-mean1) %>%
  dplyr::select(line_id, genotype_id, donor, flow_date, marker, purity, intensity)

#Import genotypes
SNPRelate::snpgdsVCF2GDS("CD14_cis_region.vcf", "CD14_cis_region1.gds",method = "copy.num.of.ref")
cd14_cis_region = gdsToMatrix("CD14_cis_region1.gds")

#CD16 (FCGR3B)
#"bcftools view -r 1:161732020 -S macrophage-gxe-study/data/sample_lists/flow_cytometry_gt_list.txt genotypes/GRCh38/imputed_20151005/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr1.GRCh38.sorted.vcf.gz > FCGR3B_qtl.vcf"
GWASTools::convertVcfGds("FCGR3B_qtl.vcf", "FCGR3B_qtl.gds")
gds <- GdsGenotypeReader("FCGR3B_qtl.gds")
cd16_qtl = data_frame(genotype_id = getVariable(gds, "sample.name"), gt = getGenotype(gds))
colnames(cd16_qtl)[colnames(cd16_qtl) == "gt"] = getVariable(gds, "snp.rs.id")
close(gds)

#Analyse CD14 only
geno_mat = as.data.frame(t(cd14_cis_region$genotypes)) %>% 
  dplyr::mutate(genotype_id = colnames(cd14_cis_region$genotypes))
cd14_data = dplyr::filter(flow_data, marker == "CD14") %>%
  dplyr::left_join(geno_mat, by = "genotype_id")
cd14_model_1 = lmer(intensity ~ (1|flow_date) + (1|donor), cd14_data)
cd14_variance = seqUtils::varianceExplained(cd14_model1)
cd14_model_2 = lmer(intensity ~ (1|flow_date) + (1|donor) + (1|rs778583), cd14_data)
cd14_variance = seqUtils::varianceExplained(cd14_model_2)

#Make a QTL plot
cd14_qtl_plot = ggplot(cd14_data, aes(x = factor(rs778583), y = intensity)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = .1))


#Analyse CD16 only
cd16_data = dplyr::filter(flow_data, marker == "CD16") %>%
  dplyr::left_join(cd16_qtl, by = "genotype_id")
cd16_model_1 = lmer(intensity ~ (1|flow_date) + (1|donor), cd16_data, REML = TRUE)
cd16_variance_1 = seqUtils::varianceExplained(cd16_model_1)
cd16_model_2 = lmer(intensity ~ (1|flow_date) + (1|donor) + (1|rs10917809), cd16_data, REML = TRUE)
cd16_variance_2 = seqUtils::varianceExplained(cd16_model_2)

#Make boxplot
cd16_qtl_plot = ggplot(cd16_data, aes(x = factor(rs10917809), y = intensity)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = .1))


#Analyse CD206 only
cd206_data = dplyr::filter(flow_data, marker == "CD206")
cd206_model = lmer(intensity ~ (1|flow_date) + (1|donor), cd206_data)
cd206_variance = seqUtils::varianceExplained(cd206_model)
