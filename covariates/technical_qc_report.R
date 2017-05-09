library("dplyr")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")

#This file makes a bunch of QC figures for the technical paper.

#Import expression and ATAC datasets
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
final_lines = dplyr::bind_rows(dplyr::select(combined_expression_data$sample_metadata, line_id, replicate), 
                               dplyr::select(atac_data$sample_metadata, line_id, replicate)) %>%
  unique()

#Set up a single output folder for all of the figures
figure_folder = "figures/supplementary/"

#Import metadata from disk
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")
line_metadata_success = dplyr::filter(line_metadata, status == "Success")

#Also import RNA RIN scores
rin1 = read.table("macrophage-gxe-study/data/bioanalyzer/bioanalyser_combined_results.txt", header = TRUE)
rin2 = read.table("macrophage-gxe-study/data/bioanalyzer/bioanalyser_illumina_bespoke.txt", header = TRUE)
rin_scores = rbind(rin1, rin2)

#How long does it take to expand iPS cells before differentiation?
duration_df = dplyr::mutate(line_metadata, received_as = ifelse(received_as == "alive/frozen", "frozen",received_as)) %>%
  dplyr::mutate(media_type = ifelse(media == "E8", "feeder-free","feeder-dependent"))
ips_duration_plot = ggplot(duration_df, aes(x = ips_culture_days, fill = received_as)) + 
  geom_histogram(binwidth = 2) + 
  facet_wrap(~media_type) +
  xlab("Duration of iPS cell culture (days)") + 
  theme_light()
ggsave(paste0(figure_folder, "diff_ips_culture_duration.pdf"), plot = ips_duration_plot, width = 7, height = 4.5)

#How long does it take from EB formation to stimulation experiment? (From EB formation to Salmonella infection)
duration_df = dplyr::mutate(line_metadata_success, differation_duration = as.numeric(salmonella_date - diff_start))
median_duration = median(duration_df$differation_duration, na.rm = TRUE)
diff_days_plot = ggplot(duration_df, aes(x = differation_duration)) + 
  geom_histogram(binwidth = 3) +
  xlab("Duration of differentiation (days)") + 
  theme_light() +
  geom_vline(xintercept = median_duration, color = "red")
ggsave(paste0(figure_folder, "diff_days_until_salmnoella_assay.pdf"), plot = diff_days_plot, width = 4.5, height = 4.5)

#How many medium changes do we do for successfull differentiations?
medium_changes_plot = ggplot(line_metadata_success, aes(x = medium_changes)) + 
  geom_histogram(binwidth = 1) + 
  xlab("Number of medium changes per differentiation") +
  theme_light()
ggsave(paste(figure_folder, "medium_changes.pdf"), plot = medium_changes_plot, width = 6, height = 5)

#Look at how many differentiations succeed and how many fail
success_df = line_metadata %>%
  dplyr::mutate(success = ifelse(status == "Success", "Success", "Fail")) %>%
  dplyr::mutate(status = ifelse(status == "RNA_QC_fail", "Stimulation_fail", status))

#What is the total success rate?
table(success_df$status)["Success"]/length(success_df$status)
#Total success is 74/117 = 0.63

#How are successes and failures distributed over time?
monthly_df = success_df %>%
  dplyr::select(donor, clone, EB_formation, status) %>% 
  dplyr::arrange(EB_formation) %>% 
  dplyr::mutate(month = format(EB_formation, "%b,%y")) %>%
  dplyr::mutate(month = factor(month, levels = unique(month)))
monthly_success_plot = ggplot(monthly_df, aes(x = month, fill = status)) + geom_bar() + 
  xlab("Start of differentiation (month)") +
  theme_light() +
  theme(axis.text.x=element_text(angle = 45))
ggsave(paste(figure_folder, "monthly_success_plot.pdf"), plot = monthly_success_plot, width = 8, height = 5)

#Look at how were the flow cytometry purity scores distributed
#Keep only lines where flow cytometry was performed within 2 weeks of salmonella infection
flow_df = dplyr::semi_join(line_metadata, final_lines, by = c("line_id", "replicate")) %>%
  dplyr::filter(!is.na(max_purity))

flow_purity_plot = ggplot(flow_df, aes(x = mean_purity)) + 
  geom_histogram(binwidth = 0.01) + 
  xlab("Mean purity") + 
  theme_light() +
  theme(legend.position = "top")
ggsave(paste0(figure_folder, "diff_flow_purity_histogram.pdf"), plot = flow_purity_plot, width = 4, height = 4)

flow_purity_plot = ggplot(flow_df, aes(x = max_purity)) + 
  geom_histogram(binwidth = 0.01) + 
  xlab("Max purity") + 
  theme_light() +
  theme(legend.position = "top")
ggsave(paste0(figure_folder, "diff_flow_purity_histogram.max.png"), plot = flow_purity_plot, width = 4, height = 4)

#Import individual flow measurements
purity_values = readRDS("macrophage-gxe-study/data/covariates/flow_cytometry_purity.rds") %>%
  dplyr::select(donor, flow_date, channel, purity) %>%
  dplyr::semi_join(flow_df, by = c("donor","flow_date"))
purity_df = tidyr::spread(purity_values, channel, purity)

cd14_cd16_scatter = ggplot(purity_df, aes(x = Pacific.Blue.A, y = PE.A)) + 
  geom_point() + 
  theme_light() + xlab("CD14 positive fraction") + 
  ylab("CD16 positive fraction")
ggsave(paste0(figure_folder, "CD14_CD16_scatter.pdf"), cd14_cd16_scatter, width = 3.5, height = 3.5)

cd14_cd206_scatter = ggplot(purity_df, aes(x = Pacific.Blue.A, y = APC.A)) + 
  geom_point() + 
  theme_light() + xlab("CD14 positive fraction") + 
  ylab("CD206 positive fraction")
ggsave(paste0(figure_folder, "CD14_CD206_scatter.pdf"), cd14_cd206_scatter, width = 3.5, height = 3.5)


#Look at sources of variation in the RNA-Seq data
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")

#How variable is RNA integrity between samples?
rin_scores_filtered = dplyr::semi_join(rin_scores, combined_expression_data$sample_metadata, by = "sample_id")
rin_plot = ggplot(rin_scores_filtered, aes(x = RIN - 0.01)) + 
  geom_histogram(binwidth = 0.2) + 
  xlab("RNA Integrity Number (RIN)") + 
  theme_light()
ggsave(paste0(figure_folder, "diff_RIN_histogram.pdf"), plot = rin_plot, width = 4.5, height = 4.5)

#iPS passage number
median_passage = median(combined_expression_data$sample_metadata$passage_diff, na.rm = TRUE)
passage_plot = ggplot(combined_expression_data$sample_metadata, aes(x = passage_diff)) + 
  geom_histogram(binwidth = 3) +
  xlab("IPSC passage number") + 
  theme_light() +
  geom_vline(xintercept = median_passage, color = "red")
ggsave(paste0(figure_folder, "diff_passage_histogram.pdf"), plot = passage_plot, width = 4.5, height = 4.5)

#iPS passage number
median_passage = median(combined_expression_data$sample_metadata$passage_diff, na.rm = TRUE)
passage_plot = ggplot(combined_expression_data$sample_metadata, aes(x = passage_diff)) + 
  geom_histogram(binwidth = 3) +
  xlab("IPSC passage number") + 
  theme_light() +
  geom_vline(xintercept = median_passage, color = "red")
ggsave(paste0(figure_folder, "diff_passage_histogram.pdf"), plot = passage_plot, width = 4.5, height = 4.5)

#iPS passage number
median_rna = median(combined_expression_data$sample_metadata$ng_ul_mean, na.rm = TRUE)
rna_plot = ggplot(combined_expression_data$sample_metadata, aes(x = ng_ul_mean)) + 
  geom_histogram(binwidth = 30) +
  xlab("RNA concentration (ng/ul)") + 
  theme_light() +
  geom_vline(xintercept = median_rna, color = "red")
ggsave(paste0(figure_folder, "diff_rna_concentration_histogram.pdf"), plot = rna_plot, width = 4.5, height = 4.5)



