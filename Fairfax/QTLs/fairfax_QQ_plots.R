library("devtools")
library("plyr")
library("ggplot2")
library("SummarizedExperiment")
library("dplyr")
load_all("../seqUtils/")
load_all("~/software/rasqual/rasqualTools/")

#Import Fairfax QTL min p-values
qtl_min_pvalues = readRDS("results/Fairfax/fairfax_qtl_min_pvalues.rds")


#Make Q-Q plots using the permutation p-values from the fairfax dataset
full_perm = purrr::map_df(qtl_min_pvalues$full, ~dplyr::mutate(.,p_eigen = p_beta, method = "Full dataset") %>% addExpectedPvalue(), identity, .id = "condition_name") %>%
  dplyr::mutate(method = "Full dataset")
shared_perm = purrr::map_df(qtl_min_pvalues$shared_84, ~dplyr::mutate(.,p_eigen = p_beta, method = "Full dataset") %>% addExpectedPvalue(), identity, .id = "condition_name") %>%
  dplyr::mutate(method = "84 samples")
qtl_df = dplyr::bind_rows(full_perm, shared_perm)

qqplot = ggplot(qtl_df, aes(x = -log(p_expected,10), y = -log(p_eigen, 10),color = method)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  facet_wrap(~condition_name, ncol = 4) +
  theme_light() + 
  ylab(expression(paste(Log[10], " observed p-value", sep = ""))) + 
  xlab(expression(paste(Log[10], " expected p-value", sep = ""))) +
  scale_colour_manual(values = c("#e66101","#5e3c99"), name = "")
ggsave("figures/supplementary/fairfax_permutation_QQ_plot.png", plot = qqplot, width = 8, height = 3.5)


#Perform eigenMT correction
eigen_MT_list = list(readr::read_tsv("processed/Fairfax/eigenMT/CD14.output.txt", col_names = TRUE),
                     readr::read_tsv("processed/Fairfax/eigenMT/IFN.output.txt", col_names = TRUE),
                     readr::read_tsv("processed/Fairfax/eigenMT/LPS2.output.txt", col_names = TRUE),
                     readr::read_tsv("processed/Fairfax/eigenMT/LPS24.output.txt", col_names = TRUE)) %>%
  purrr::map_df(identity)
eigenMT_results = eigen_MT_list %>%
  dplyr::transmute(phenotype_id = gene, n_eigen_snps = TESTS) %>%
  dplyr::filter(phenotype_id != "gene") %>%
  dplyr::mutate(n_eigen_snps = as.integer(n_eigen_snps)) %>%
  dplyr::distinct()

full_eigen = dplyr::left_join(full_perm, eigenMT_results, by = "phenotype_id") %>% dplyr::group_by(phenotype_id, condition_name) %>%
  dplyr::filter(!is.na(n_eigen_snps)) %>%
  dplyr::mutate(p_eigen = p.adjust(p_nominal, method = "bonferroni", n = n_eigen_snps)) %>%
  dplyr::group_by(condition_name) %>%
  dplyr::mutate(p_eigen_fdr = p.adjust(p_eigen, "fdr")) %>%
  addExpectedPvalue(.) %>%
  ungroup()
shared_eigen = dplyr::left_join(shared_perm, eigenMT_results, by = "phenotype_id") %>% dplyr::group_by(phenotype_id, condition_name) %>%
  dplyr::filter(!is.na(n_eigen_snps)) %>%
  dplyr::mutate(p_eigen = p.adjust(p_nominal, method = "bonferroni", n = n_eigen_snps)) %>%
  dplyr::group_by(condition_name) %>%
  dplyr::mutate(p_eigen_fdr = p.adjust(p_eigen, "fdr")) %>%
  addExpectedPvalue(.) %>%
  ungroup()
qtl_df = dplyr::bind_rows(full_eigen, shared_eigen)

qqplot = ggplot(qtl_df, aes(x = -log(p_expected,10), y = -log(p_eigen, 10),color = method)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  facet_wrap(~condition_name, ncol = 4) +
  theme_light() + 
  ylab(expression(paste(Log[10], " observed p-value", sep = ""))) + 
  xlab(expression(paste(Log[10], " expected p-value", sep = ""))) +
  scale_colour_manual(values = c("#e66101","#5e3c99"), name = "")
ggsave("figures/supplementary/fairfax_eigenMT_QQ_plot.png", plot = qqplot, width = 8, height = 3.5)



#Perform eigenMT correction
full_perm = purrr::map_df(qtl_min_pvalues$full, ~dplyr::mutate(.,p_eigen = p_beta, method = "Full dataset") %>% addExpectedPvalue(), identity, .id = "condition_name") %>%
  dplyr::mutate(method = "Full dataset")
shared_perm = purrr::map_df(qtl_min_pvalues$shared_42, ~dplyr::mutate(.,p_eigen = p_beta, method = "Full dataset") %>% addExpectedPvalue(), identity, .id = "condition_name") %>%
  dplyr::mutate(method = "42 samples")

eigen_MT_list = list(readr::read_tsv("processed/Fairfax/eigenMT/CD14.output.txt", col_names = TRUE),
                     readr::read_tsv("processed/Fairfax/eigenMT/IFN.output.txt", col_names = TRUE),
                     readr::read_tsv("processed/Fairfax/eigenMT/LPS2.output.txt", col_names = TRUE),
                     readr::read_tsv("processed/Fairfax/eigenMT/LPS24.output.txt", col_names = TRUE)) %>%
  purrr::map_df(identity)
eigenMT_results = eigen_MT_list %>%
  dplyr::transmute(phenotype_id = gene, n_eigen_snps = TESTS) %>%
  dplyr::filter(phenotype_id != "gene") %>%
  dplyr::mutate(n_eigen_snps = as.integer(n_eigen_snps)) %>%
  dplyr::distinct()

full_eigen = dplyr::left_join(full_perm, eigenMT_results, by = "phenotype_id") %>% dplyr::group_by(phenotype_id, condition_name) %>%
  dplyr::filter(!is.na(n_eigen_snps)) %>%
  dplyr::mutate(p_eigen = p.adjust(p_nominal, method = "bonferroni", n = n_eigen_snps)) %>%
  dplyr::group_by(condition_name) %>%
  dplyr::mutate(p_eigen_fdr = p.adjust(p_eigen, "fdr")) %>%
  addExpectedPvalue(.) %>%
  ungroup()
shared_eigen = dplyr::left_join(shared_perm, eigenMT_results, by = "phenotype_id") %>% dplyr::group_by(phenotype_id, condition_name) %>%
  dplyr::filter(!is.na(n_eigen_snps)) %>%
  dplyr::mutate(p_eigen = p.adjust(p_nominal, method = "bonferroni", n = n_eigen_snps)) %>%
  dplyr::group_by(condition_name) %>%
  dplyr::mutate(p_eigen_fdr = p.adjust(p_eigen, "fdr")) %>%
  addExpectedPvalue(.) %>%
  ungroup()
qtl_df = dplyr::bind_rows(full_eigen, shared_eigen)

qqplot = ggplot(qtl_df, aes(x = -log(p_expected,10), y = -log(p_eigen, 10),color = method)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  facet_wrap(~condition_name, ncol = 4) +
  theme_light() + 
  ylab(expression(paste(Log[10], " observed p-value", sep = ""))) + 
  xlab(expression(paste(Log[10], " expected p-value", sep = ""))) +
  scale_colour_manual(values = c("#e66101","#5e3c99"), name = "")
ggsave("figures/supplementary/fairfax_eigenMT_QQ_plot_n42.png", plot = qqplot, width = 8, height = 3.5)