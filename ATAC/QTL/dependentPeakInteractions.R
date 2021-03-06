library("purrr")
library("ggplot2")
library("devtools")
library("dplyr")
library("wiggleplotr")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")
load_all("macrophage-gxe-study/housekeeping/")

#Import ATAC data
atac_data = readRDS("results/ATAC/ATAC_combined_accessibility_data_covariates.rds")

#Import the VCF file
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")

#List of ATAC tabix files
atac_tabix_list = qtlResults()$atac_rasqual

#Import minimal p-values
rasqual_min_pvalues = readRDS("results/ATAC/QTLs/rasqual_min_pvalues.rds")
min_pvalue_hits = lapply(rasqual_min_pvalues, function(x){dplyr::filter(x, p_eigen < fdr_thresh)})
min_pvalues_df = purrr::map_df(min_pvalue_hits, identity, .id = "condition_name") %>%
  dplyr::arrange(gene_id, p_nominal)

#Extract data
naive_ifng_atac = extractConditionFromExpressionList(c("naive","IFNg"), atac_data)
naive_sl1344_atac = extractConditionFromExpressionList(c("naive","SL1344"), atac_data)
naive_ifng_sl1344_atac = extractConditionFromExpressionList(c("naive","IFNg_SL1344"), atac_data)

#Specify models of interest
model1 = as.formula("cqn ~genotype + peak_type + condition_name + cqn_PC1 + cqn_PC2 + cqn_PC3 + peak_type*condition_name + 
               genotype*peak_type + genotype*condition_name + genotype*condition_name*peak_type")
model0 = as.formula("cqn ~genotype + peak_type + condition_name  + cqn_PC1 + cqn_PC2 + cqn_PC3 + peak_type*condition_name + 
                    genotype*peak_type + genotype*condition_name")

#Import list of dependent peaks
result_list = readRDS("results/ATAC/QTLs/qtl_peak_type_assignment.rds")
dependent_peaks = result_list$dependents

#Filter the VCF file for quicker testing
genotypes = vcf_file$genotypes[unique(dependent_peaks$unique_masters$snp_id),]
snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)

#Test peak-peak-genotype interactions in naive_vs_IFNg setting
ifng_interactions = purrr::by_row(dependent_peaks$unique_masters, testThreewayInteraction, naive_ifng_atac$cqn, 
                  naive_ifng_atac$sample_metadata, vcf_file, model0, model1, .collate = "rows", .to = "p_nominal")
sl1344_interactions = purrr::by_row(dependent_peaks$unique_masters, testThreewayInteraction, naive_sl1344_atac$cqn, 
                                    naive_sl1344_atac$sample_metadata, vcf_file, model0, model1, .collate = "rows", .to = "p_nominal")
ifng_sl1344_interactions = purrr::by_row(dependent_peaks$unique_masters, testThreewayInteraction, naive_ifng_sl1344_atac$cqn, 
                                         naive_ifng_sl1344_atac$sample_metadata, vcf_file, model0, model1, .collate = "rows", .to = "p_nominal")
interaction_list = list(IFNg = ifng_interactions, SL1344 = sl1344_interactions, IFNg_SL1344 = ifng_sl1344_interactions)
saveRDS(interaction_list, "results/ATAC/QTLs/peak_peak_interactions.txt")
interaction_list = readRDS("results/ATAC/QTLs/peak_peak_interactions.txt")

#Permute genotypes
new_labels = sample(colnames(vcf_file$genotypes))
vcf_file_perm = vcf_file
colnames(vcf_file_perm$genotypes) = new_labels

#Test interactions with permuted genotypes
ifng_interactions = purrr::by_row(dependent_peaks$unique_masters, testThreewayInteraction, naive_ifng_atac$cqn, 
                                  naive_ifng_atac$sample_metadata, vcf_file_perm, model0, model1, .collate = "rows", .to = "p_nominal")
sl1344_interactions = purrr::by_row(dependent_peaks$unique_masters, testThreewayInteraction, naive_sl1344_atac$cqn, 
                                    naive_sl1344_atac$sample_metadata, vcf_file_perm, model0, model1, .collate = "rows", .to = "p_nominal")
ifng_sl1344_interactions = purrr::by_row(dependent_peaks$unique_masters, testThreewayInteraction, naive_ifng_sl1344_atac$cqn, 
                                         naive_ifng_sl1344_atac$sample_metadata, vcf_file_perm, model0, model1, .collate = "rows", .to = "p_nominal")
interaction_list = list(IFNg = ifng_interactions, SL1344 = sl1344_interactions, IFNg_SL1344 = ifng_sl1344_interactions)
saveRDS(interaction_list, "results/ATAC/QTLs/peak_peak_interactions.permuted.txt")
interaction_list_perm = readRDS("results/ATAC/QTLs/peak_peak_interactions.permuted.txt")

#Alternative option - permute conditions withing individual (less bad, but still does not work):
#Permute conditions
perm_conditions = dplyr::group_by(naive_ifng_atac$sample_metadata, donor) %>% 
  dplyr::mutate(condition_char_new = sample(condition_char)) %>% 
  dplyr::select(donor, condition_char, condition_char_new) %>% dplyr::ungroup() %>% 
  dplyr::mutate(perm_sample_id = paste(donor, condition_char_new, "ATAC", sep = "_"))
colnames(naive_ifng_atac$cqn) = perm_conditions$perm_sample_id

ifng_interactions = purrr::by_row(dependent_peaks$unique_masters, testThreewayInteraction, naive_ifng_atac$cqn, 
                                  naive_ifng_atac$sample_metadata, filtered_vcf, model0, model1, .collate = "rows", .to = "p_nominal")



#Compare p-value distributions for nominal and permutation runs
interaction_df = purrr::map_df(interaction_list, identity, .id = "other_condition") %>%
  dplyr::mutate(test = "nominal")
interaction_df_perm = purrr::map_df(interaction_list_perm, identity, .id = "other_condition")  %>%
  dplyr::mutate(test = "permuted")
df_all = dplyr::bind_rows(interaction_df, interaction_df_perm)

p_histogram = ggplot(interaction_df, aes(x = p_nominal)) + geom_histogram(bins = 20) + 
  #facet_wrap(~test) + 
  theme_light() +
  xlab("p-value")
ggsave("figures/supplementary/dependent_peak_interaction_histogram.pdf", plot = p_histogram, width = 5, height = 3)
ggsave("figures/supplementary/dependent_peak_interaction_histogram.png", plot = p_histogram, width = 5, height = 3)


#Convert list into a df and extract interaction hits
interaction_df = purrr::map_df(interaction_list, identity, .id = "other_condition") %>%
  dplyr::mutate(p_nominal = pmin(p_nominal*3,1)) %>%
  dplyr::group_by(master_id, dependent_id) %>%
  dplyr::arrange(p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
interaction_hits = dplyr::mutate(interaction_df, p_fdr = p.adjust(p_nominal, method = "fdr")) %>% 
  dplyr::arrange(p_nominal) %>% dplyr::filter(p_fdr < 0.1) %>%
  dplyr::mutate(baseline_condition = "naive") 

#Make a Q-Q plot for the interaction test
qq_df = dplyr::mutate(interaction_df, p_eigen = p_nominal) %>% addExpectedPvalue()
qq_plot = ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_nominal,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")

#Extract RASQUAL effect sizes for all peak-peak pairs
unique_snps = interaction_hits$snp_id
snp_ranges = dplyr::filter(vcf_file$snpspos, snpid %in% unique_snps) %>% 
  dplyr::transmute(snp_id = snpid, seqnames = chr, start = pos, end = pos, strand = "+") %>% 
  dataFrameToGRanges()
snp_results = purrr::map_df(atac_tabix_list, ~rasqualTools::tabixFetchSNPs(snp_ranges, .),.id = "condition_name")

#Extract effect sizes for each gene_snp pair
inter_hits = dplyr::select(interaction_hits, master_id, dependent_id, snp_id, p_nominal, baseline_condition, other_condition)
snp_effects = dplyr::select(snp_results, gene_id, snp_id, condition_name, beta)

effects = dplyr::left_join(inter_hits, snp_effects, by = c("snp_id","master_id" = "gene_id", "baseline_condition" = "condition_name")) %>% 
  dplyr::rename(master_baseline = beta) %>%
  dplyr::left_join(snp_effects, by = c("snp_id","master_id" = "gene_id", "other_condition" = "condition_name")) %>% 
  dplyr::rename(master_other = beta)  %>%
  dplyr::left_join(snp_effects, by = c("snp_id","dependent_id" = "gene_id", "baseline_condition" = "condition_name")) %>% 
  dplyr::rename(dependent_baseline = beta) %>% 
  dplyr::left_join(snp_effects, by = c("snp_id","dependent_id" = "gene_id", "other_condition" = "condition_name")) %>% 
  dplyr::rename(dependent_other = beta)

#Calculate diffs
effects_diff = dplyr::mutate(effects, master_diff = master_other - master_baseline, dependent_diff = dependent_other - dependent_baseline)

#Filter interaction results by effect size
filtered_interactions = dplyr::filter(effects_diff, abs(master_baseline) > 0.59, abs(dependent_diff) > 0.59, abs(master_diff) < 1) %>% 
  dplyr::group_by(master_id, dependent_id) %>% 
  arrange(p_nominal) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()

#Extract raw data for the filtered events
ifng_effects = purrr::by_row(dplyr::filter(filtered_interactions, other_condition == "IFNg"), testThreewayInteraction, naive_ifng_atac$cqn, 
                                  naive_ifng_atac$sample_metadata, vcf_file, model0, model1, p_only = FALSE)
sl1344_effects = purrr::by_row(dplyr::filter(filtered_interactions, other_condition == "SL1344"), testThreewayInteraction, naive_sl1344_atac$cqn, 
                                    naive_sl1344_atac$sample_metadata, vcf_file, model0, model1, p_only = FALSE)
ifng_sl1344_effects = purrr::by_row(dplyr::filter(filtered_interactions, other_condition == "IFNg_SL1344"), testThreewayInteraction, naive_ifng_sl1344_atac$cqn, 
                                         naive_ifng_sl1344_atac$sample_metadata, vcf_file, model0, model1, p_only = FALSE)
#Make plots for all of the hits
joint_data = bind_rows(ifng_effects, sl1344_effects, ifng_sl1344_effects)
joint_list = purrr::map(as.list(joint_data$.out),function(x){x$data})
plot_list = map(joint_list, ~ggplot(., aes(x = factor(genotype), y = cqn)) + geom_point() + geom_boxplot() + facet_grid(peak_type ~ condition_name))
names(plot_list) = paste(joint_data$dependent_id, joint_data$master_id, sep = "-")
savePlotList(plot_list, "results/ATAC/QTLs/peak-peak_interactions/selected_interactions/", suffix = ".pdf", width = 8, height = 8)

##### Distances between master and dependent peaks ####
#Calculated distances between master and dependent peaks
peak_distances = dplyr::select(dependent_peaks$unique_masters, dependent_id, master_id) %>%
  calculatePeakDistance(atac_data$gene_metadata) %>%
  dplyr::select(dependent_id, master_id, distance)
peak_distances_cl = dplyr::select(dependent_peaks$cluster_master, dependent_id, master_id) %>%
  calculatePeakDistance(atac_data$gene_metadata) %>%
  dplyr::select(dependent_id, master_id, distance)

#Make a histogram of distances
peak_distance_df = rbind(peak_distances, peak_distances_cl)
plot = ggplot(peak_distance_df, aes(x = distance/1000)) + geom_histogram(binwidth = 1) + 
  xlab("Distance between master and dependent peaks (kb)") + 
  ylab("Number of peak pairs")
ggsave("results/ATAC/QTLs/properties/master_dependent_peak_distance.pdf", plot = plot, width = 6, height = 6)


#Calculate distances within peak clusters
peak_clusters = readRDS("results/ATAC/QTLs/clustered_qtl_peaks.rds")

#Construct master-dependent pairs
dependents = dplyr::transmute(peak_clusters$cluster_memberships, dependent_id = gene_id, cluster_id)
masters = dplyr::transmute(peak_clusters$lead_snps, master_id = gene_id, cluster_id)
cluster_distances = dplyr::left_join(dependents, masters, by = "cluster_id") %>% 
  dplyr::filter(master_id != dependent_id) %>% calculatePeakDistance(atac_data$gene_metadata)
plot = ggplot(dplyr::filter(cluster_distances, abs(distance) < 50000), aes(x = distance/1000)) + geom_histogram(binwidth = 1) + 
  xlab("Distance between cluster lead peak and dependent peaks (kb)") + 
  ylab("Number of peak pairs")
ggsave("results/ATAC/QTLs/properties/cluster_dependent_peak_distance.pdf", plot = plot, width = 6, height = 6)



#Make coverage plots of master and dependent peak pairs

#Construct metadata df for wiggleplotr
meta = wiggleplotrConstructMetadata(atac_data$counts, atac_data$sample_metadata, "/Volumes/Ajamasin/bigwig/ATAC/")

#Add peak coords to the dependent pairs
peak_coords = dplyr::select(atac_data$gene_metadata, gene_id, chr, start, end) %>% tbl_df()
joint_coords = dplyr::left_join(joint_data, peak_coords, by = c("master_id" = "gene_id")) %>%
  dplyr::left_join(peak_coords, by = c("dependent_id" = "gene_id"))

#Construct df regions
joint_regions = dplyr::transmute(joint_coords, master_id, dependent_id, snp_id, p_nominal, chr = chr.x, 
                                 region_start = pmin(start.x, start.y), region_end = pmax(end.x, end.y),
                                 max_fc = pmax(master_baseline,master_other,dependent_baseline,dependent_other),
                                 min_fc = pmin(master_baseline,master_other,dependent_baseline,dependent_other)) %>%
  dplyr::mutate(max_sign = ifelse(abs(max_fc) > abs(min_fc), max_fc, min_fc)) %>%
  dplyr::arrange(p_nominal)


coverageByRow <- function(row_df, meta_df, gene_metadata, genotypes){
  
  #Extract peaks from row
  peaks_df = dplyr::filter(gene_metadata, chr == row_df$chr, start >= row_df$region_start, end <= row_df$region_end)
  peak_annot = wiggpleplotrConstructPeakAnnotations(peaks_df)
  track_data = wiggleplotrGenotypeColourGroup(meta_df, row_df$snp_id, genotypes, row_df$max_sign)
  print(row_df$master_id)
  
  #Make coverage plot
  coverage = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = track_data, rescale_introns = FALSE, 
                          transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(), flanking_length = c(500,500), 
                          connect_exons = FALSE, label_type = "peak", plot_fraction = 0.2, heights = c(0.8,0.2), coverage_type = "line")
  return(coverage)
}

#Make a lot of coverage plots
coverage_df = purrr::by_row(joint_regions, ~coverageByRow(.,meta, atac_data$gene_metadata, vcf_file$genotypes))
plot_list = coverage_df$.out
names(plot_list) = paste(coverage_df$master_id, coverage_df$dependent_id, sep = "-")

#Save them to disk
savePlotList(plot_list, "results/ATAC/QTLs/peak-peak_interactions/selected_coverage/", width = 5, height = 6, suffix = ".png")

#Make single example plot
coverage_df = purrr::by_row(joint_regions[1,], ~coverageByRow(.,meta, atac_data$gene_metadata, vcf_file$genotypes))

