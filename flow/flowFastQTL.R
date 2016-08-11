library("devtools")
load_all("../seqUtils/")
library("lme4")


#Import data
flow_purity = readRDS("macrophage-gxe-study/data/covariates/flow_cytometry_purity.rds")
line_medatada = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")

#Extract gene coordinates from gene expression data
eqtl_data_list = readRDS("results/SL1344/combined_expression_data_covariates.rds")
gene_coords = dplyr::filter(eqtl_data_list$gene_metadata, gene_id %in% 
                              c("ENSG00000170458","ENSG00000203747","ENSG00000162747","ENSG00000260314")) %>%
  dplyr::select(-gene_id) %>%
  dplyr::rename(gene_id = gene_name) %>%
  dplyr::mutate(gene_name = gene_id)

#Import genotypes
SNPRelate::snpgdsVCF2GDS("flow/genotypes/flow_cis_regions.vcf", "flow/genotypes/flow_cis_regions.gds",method = "copy.num.of.ref")
flow_cis_region = gdsToMatrix("flow/genotypes/flow_cis_regions.gds")

#Map channel to marker
channel_marker_map = data_frame(channel = c("APC.A","PE.A","PE.A","Pacific.Blue.A"), 
                                gene_name = c("MRC1","FCGR3A","FCGR3B","CD14"))
unique_lines = dplyr::select(line_medatada, line_id, donor, genotype_id) %>% unique()
flow_data = dplyr::left_join(flow_purity, channel_marker_map, by = "channel") %>%
  dplyr::mutate(donor = ifelse(donor == "fpdj", "nibo",donor)) %>% #fpdj and nibo are the same donors
  dplyr::left_join(unique_lines, by = "donor") %>%
  dplyr::mutate(intensity = mean2-mean1) %>%
  dplyr::select(line_id, genotype_id, donor, flow_date, gene_name, purity, intensity)

#Construct a matrix of intensity values
flow_df = dplyr::select(flow_data,line_id, genotype_id, flow_date, gene_name, intensity) %>% 
  tidyr::spread(gene_name, intensity) %>%
  dplyr::mutate(sample_id = paste(flow_df$line_id, as.character(flow_df$flow_date), sep = "_"))

#Make a matrix of flow data and perform PCA
flow_matrix = as.matrix(flow_df[,c(4,5,7)])
rownames(flow_matrix) = flow_df$sample_id
pca_res = prcomp(flow_matrix, scale = TRUE, centre = TRUE)

#Choose outliers based on PCA and remove them
outlier_samples = c("hayt_1_2015-10-16","fafq_1_2015-10-16","iill_1_2015-10-20")
flow_df_filtered = dplyr::filter(flow_df, !(sample_id %in% outlier_samples))

#Calculate mean intensities for QTL mapping
flow_mean_df = dplyr::group_by(flow_df_filtered, genotype_id) %>% 
  dplyr::summarise(CD14 = mean(CD14), FCGR3A = mean(FCGR3A), FCGR3B = mean(FCGR3B), MRC1 = mean(MRC1))
flow_mean_matrix = t(as.matrix(flow_mean_df[,2:5]))
colnames(flow_mean_matrix) = flow_mean_df$genotype_id

#Random sample per donor
flow_random_df = dplyr::group_by(flow_df_filtered, genotype_id) %>% dplyr::filter(row_number() == 1)
flow_random_matrix = t(as.matrix(flow_random_df[,4:7]))
colnames(flow_random_matrix) = flow_random_df$genotype_id

#Export intensity values for FastQTL
fastqtl_gene_pos = constructFastQTLGenePos(gene_coords)
mean_exp = prepareFastqtlMatrix(flow_mean_matrix, fastqtl_gene_pos)
random_exp = prepareFastqtlMatrix(flow_random_matrix, fastqtl_gene_pos)
write.table(mean_exp, "results/flow/fastqtl/mean_intensity.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(random_exp, "results/flow/fastqtl/random_intensity.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

#Calculate principal components
pca_random = prcomp(t(flow_random_matrix[c(1,2,4),]), scale = TRUE, centre = TRUE)
pca_mean = prcomp(t(flow_mean_matrix[c(1,2,4),]), scale = TRUE, centre = TRUE)

#Make cov matrices
random_cov = t(pca_random$x[,1, drop = F]) %>% as.data.frame() %>% 
  dplyr::mutate(id = "PC1") %>% dplyr::select(id, everything())
write.table(random_cov, "results/flow/fastqtl/random_cov.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

mean_cov = t(pca_mean$x[,1, drop = F]) %>% as.data.frame() %>% 
  dplyr::mutate(id = "PC1") %>% dplyr::select(id, everything())
write.table(mean_cov, "results/flow/fastqtl/mean_cov.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

#Make some QTL plots
flow_random_df = dplyr::mutate(as.data.frame(t(flow_random_matrix)), genotype_id = colnames(flow_random_matrix))
genotype_df = dplyr::mutate(as.data.frame(t(flow_cis_region$genotypes)), genotype_id = colnames(flow_cis_region$genotypes))
plot_df = dplyr::left_join(flow_random_df, genotype_df, by = "genotype_id") %>%
  dplyr::mutate(rs4657019 = as.factor(rs4657019)) %>%
  dplyr::mutate(rs2569177 = as.factor(rs2569177))

cd16_plot = ggplot(plot_df, aes(x = rs4657019, y = FCGR3A)) + 
  geom_boxplot() +
  ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .1)) + 
  theme_light() + 
  ylab("CD16 intensity")
ggsave("figures/supplementary/flow_CD16_QTL.pdf", plot = cd16_plot, width = 3.5, height = 3.5)

cd14_plot = ggplot(plot_df, aes(x = rs2569177, y = CD14)) + 
  geom_boxplot() +
  ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .1)) + 
  theme_light() + 
  ylab("CD14 intensity")
ggsave("figures/supplementary/flow_CD14_QTL.pdf", plot = cd14_plot, width = 3.5, height = 3.5)


#Quick variance component analysis
#Estimate variance explained by different factors
cd14_exp = lmer(CD14 ~ (1|flow_date) + (1|line_id), flow_df_filtered) %>% varianceExplained()
cd16_exp = lmer(FCGR3A ~ (1|flow_date) + (1|line_id), flow_df_filtered) %>% varianceExplained()
cd206_exp = lmer(MRC1 ~ (1|flow_date) + (1|line_id), flow_df_filtered) %>% varianceExplained()

#Estimate variance explained
marker_var_exp = dplyr::bind_rows(cd14_exp, cd16_exp, cd206_exp) %>% 
  dplyr::mutate(marker = c("CD14","CD16","CD206")) %>% 
  dplyr::select(-type) %>%
  dplyr::rename(date = flow_date, line = line_id, residual = Residual) %>%
  tidyr::gather(source, variance, date:residual) %>%
  dplyr::mutate(source = factor(source, levels = c("residual", "date", "line", "QTL"))) %>%
  dplyr::arrange(source)

#Make a plot of variance explained
var_explained_1 = ggplot(marker_var_exp, aes(x = marker, y = variance, fill = source)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  xlab("Surface marker") + 
  ylab("Proportion variance explained") + 
  scale_fill_brewer(palette = "Set1")
ggsave("figures/supplementary/flow_var_explained.pdf", plot = var_explained_1, width = 4, height = 3.5)


#Restimate variance explaiend with the QTLs included in the model
qtl_variants = plot_df[,c("genotype_id","rs2569177","rs4657019")]
flow_df_qtl = dplyr::left_join(flow_df_filtered, qtl_variants)

#Estimate variance explained by different factors
cd14_exp = lmer(CD14 ~ (1|flow_date) + (1|line_id) + (1|rs2569177) + (1|rs4657019), flow_df_qtl) %>% varianceExplained()
cd16_exp = lmer(FCGR3A ~ (1|flow_date) + (1|line_id) + (1|rs2569177) + (1|rs4657019), flow_df_qtl) %>% varianceExplained()
cd206_exp = lmer(MRC1 ~ (1|flow_date) + (1|line_id) + (1|rs2569177) + (1|rs4657019), flow_df_qtl) %>% varianceExplained()

#Estimate variance explained
marker_var_exp2 = dplyr::bind_rows(cd14_exp, cd16_exp, cd206_exp) %>% 
  dplyr::mutate(marker = c("CD14","CD16","CD206")) %>% 
  dplyr::select(-type) %>%
  dplyr::rename(date = flow_date, line = line_id, residual = Residual) %>%
  tidyr::gather(source, variance, date:rs4657019) %>%
  dplyr::mutate(source = as.character(source)) %>%
  dplyr::mutate(source = ifelse(source %in% c("rs2569177","rs4657019"), "QTL",source)) %>%
  dplyr::mutate(source = factor(source, levels = c("residual", "date", "line", "QTL"))) %>%
  dplyr::arrange(source)

#Make a new plot with QTLs included
var_explained_2 = ggplot(marker_var_exp2, aes(x = marker, y = variance, fill = source)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  xlab("Surface marker") + 
  ylab("Proportion variance explained") +
  scale_fill_brewer(palette = "Set1")
ggsave("figures/supplementary/flow_var_explained_QTL.pdf", plot = var_explained_2, width = 4, height = 3.5)


#Some summary stats
#Extract the donors that have more than one replicate
replicate_donors = dplyr::filter(flow_df) %>% 
  group_by(line_id) %>% 
  dplyr::summarise(count = length(line_id)) %>% 
  dplyr::filter(count > 1)
