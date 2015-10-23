library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("ggplot2")

#This file makes a bunch of QC figures for the technical paper.

#Set up a single output folder for all of the figures
figure_folder = "results/technical_report/"

#Import metadata from disk
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds")
line_metadata_success = dplyr::filter(line_metadata, status == "Success")

#Also import RNA RIN scores
rin1 = read.table("macrophage-gxe-study/data/bioanalyzer/bioanalyser_combined_results.txt", header = TRUE)
rin2 = read.table("macrophage-gxe-study/data/bioanalyzer/bioanalyser_illumina_bespoke.txt", header = TRUE)
rin_scores = rbind(rin1, rin2)

#How long does it take to expand iPS cells before differentiation?
ips_duration_plot = ggplot(line_metadata, aes(x = ips_culture_days, fill = received_as)) + 
  geom_histogram(binwidth = 2) + 
  facet_wrap(~media) +
  xlab("Duration of iPS cell culture (days)")
ggsave(paste(figure_folder, "ips_culture_duration.pdf"), plot = ips_duration_plot, width = 8, height = 5)

#How long does it take from EB formation to stimulation experiment? (From EB formation to Salmonella infection)
diff_days_plot = ggplot(line_metadata_success, aes(x = diff_days)) + 
  geom_histogram(binwidth = 2) +
  xlab("Duration of differentiation (days)")
ggsave(paste(figure_folder, "salmonella_diff_duration.pdf"), plot = diff_days_plot, width = 6, height = 5)

#How many medium changes do we do for successfull differentiations?
medium_changes_plot = ggplot(line_metadata_success, aes(x = medium_changes)) + 
  geom_histogram(binwidth = 1) + 
  xlab("Number of medium changes per differentiation")
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
  theme(axis.text.x=element_text(angle = 45))
ggsave(paste(figure_folder, "monthly_success_plot.pdf"), plot = monthly_success_plot, width = 8, height = 5)

#Look at how were the flow cytometry purity scores distributed
#Keep only lines where flow cytometry was performed within 2 weeks of salmonella infection
flow_df = line_metadata %>% 
  dplyr::filter(donor != "mijn") %>%
  dplyr::filter(status %in% c("Success", "FC_QC_fail"), !is.na(max_purity)) %>%
  dplyr::filter((status == "Success" & abs(flow_date - salmonella) < 14) | status == "FC_QC_fail") 

flow_purity_plot = ggplot(flow_df, aes(x = max_purity-0.001, fill = status)) + 
  geom_histogram(binwidth = 0.01) + 
  xlab("Maximum purity")
ggsave(paste(figure_folder, "flow_purity_histogram.pdf"), plot = flow_purity_plot, width = 7, height = 5)

#How variable is the mean RNA concentration between different lines?
mean_rna_concentration = ggplot(line_metadata_success, aes(x = ng_ul_mean)) + 
  geom_histogram(binwidth = 20) + 
  xlab("Mean RNA concentration (ng/ul)")
ggsave(paste(figure_folder, "mean_rna_concentration.pdf"), plot = mean_rna_concentration, width = 6, height = 5)

#How variable is RNA integrity between samples?
rin_plot = ggplot(rin_scores, aes(x = RIN-0.001)) + 
  geom_histogram(binwidth = 0.2) + 
  xlab("RNA Integrity Number (RIN)") + 
  scale_x_continuous(limits = c(0,10))
ggsave(paste(figure_folder, "RIN_histogram.pdf"), plot = rin_plot, width = 6, height = 5)


