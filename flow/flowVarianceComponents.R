library("dplyr")
library("lme4")
library("devtools")
load_all("../seqUtils/")

#Import data
flow_purity = readRDS("macrophage-gxe-study/data/covariates/flow_cytometry_purity.rds")
line_medatada = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")
vcf_file = readRDS("genotypes/SL1344/array_genotypes.59_samples.imputed.vcfToMatrix.rds")

#Map channel to marker
channel_marker_map = data_frame(channel = c("APC.A","PE.A","Pacific.Blue.A"), marker = c("CD206","CD16","CD14"))
unique_lines = dplyr::select(line_medatada, line_id, donor, genotype_id) %>% unique()
flow_data = dplyr::left_join(flow_purity, channel_marker_map, by = "channel") %>%
  dplyr::left_join(unique_lines, by = "donor") %>%
  dplyr::mutate(donor = ifelse(donor == "fpdj", "nibo",donor)) %>% #fpdj and nibo are the same donors
  dplyr::mutate(intensity = mean2-mean1) %>%
  dplyr::select(line_id, genotype_id, donor, flow_date, marker, purity, intensity)

#Import genotypes
a = vcfToMatrix("CD14_qtl.vcf", "GRCh38", genotype_sep = "/")
cd14_qtl = dplyr::mutate(as.data.frame(t(a$genotypes)), genotype_id = colnames(a$genotypes))

#Analyse CD14 only
cd14_data = dplyr::filter(flow_data, marker == "CD14") %>%
  dplyr::left_join(cd14_qtl, by = "genotype_id")
cd14_model = lmer(intensity ~ (1|flow_date) + (1|donor), cd14_data)
cd14_variance = seqUtils::varianceExplained(cd14_model)
cd14_model = lmer(intensity ~ (1|flow_date) + (1|donor) + (1|rs778583), cd14_data)
cd14_variance = seqUtils::varianceExplained(cd14_model)


#Analyse CD16 only
cd16_data = dplyr::filter(flow_data, marker == "CD16")
cd16_model = lmer(intensity ~ (1|flow_date) + (1|donor), cd16_data)
cd16_variance = seqUtils::varianceExplained(cd16_model)

#Analyse CD206 only
cd206_data = dplyr::filter(flow_data, marker == "CD206")
cd206_model = lmer(intensity ~ (1|flow_date) + (1|donor), cd206_data)
cd206_variance = seqUtils::varianceExplained(cd206_model)