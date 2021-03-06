library("purrr")
library("gProfileR")
library("dplyr")

#Imoprt Mfuzz clustering results
clusters = readRDS("results/SL1344/DE/mFuzz_clusters_for_GO_analysis.rds")

#Run gProfilers
cluster1 = (dplyr::filter(clusters, new_cluster_id == 1) %>% dplyr::arrange(-IFNg_SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE, exclude_iea = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster2 = (dplyr::filter(clusters, new_cluster_id == 2) %>% dplyr::arrange(-SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE, exclude_iea = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster3 = (dplyr::filter(clusters, new_cluster_id == 3) %>% dplyr::arrange(-IFNg_SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE, exclude_iea = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster4 = (dplyr::filter(clusters, new_cluster_id == 4) %>% dplyr::arrange(-SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE, exclude_iea = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster5 = (dplyr::filter(clusters, new_cluster_id == 5) %>% dplyr::arrange(-IFNg_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE, exclude_iea = TRUE, significant = FALSE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster6 = (dplyr::filter(clusters, new_cluster_id == 6) %>% dplyr::arrange(IFNg_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE, exclude_iea = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

cluster789 = (dplyr::filter(clusters, new_cluster_id %in% c(7,8,9)) %>% dplyr::arrange(IFNg_SL1344_log2FC))$gene_id %>% 
  gprofiler(max_set_size = 3000, ordered_query = TRUE, exclude_iea = TRUE) %>%
  dplyr::select(term.size, query.size, overlap.size, p.value, term.name, domain, term.id) %>% 
  dplyr::arrange(p.value) %>% tbl_df()

gprofiler_results = list(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6, cluster789)
names(gprofiler_results) = c("1","2","3","4","5","6","7,8,9")
saveRDS(gprofiler_results, "results/SL1344/DE/Mfuzz_clustering_go_enrichments.rds")


#Read GO enrichments form disk
go_enrichments = readRDS("results/SL1344/DE/Mfuzz_clustering_go_enrichments.rds")
go_df = purrr::map_df(go_enrichments, identity ,.id = "cluster_id") %>%
  dplyr::mutate(cluster_id = ifelse(cluster_id == "7,8,9", "7-9", cluster_id))

selected_terms = c("Antigen processing and presentation","cell death",
                   "MAPK signaling pathway","NF-kappa B signaling pathway","Cell cycle",
                   "response to lipopolysaccharide","TNF signaling pathway","cell migration",
                   "type I interferon receptor binding","JAK-STAT cascade", "Interferon alpha/beta signaling",
                   "response to interferon-gamma","mitotic cell cycle","ncRNA processing")


selected_terms = c("Antigen processing and presentation","cell death",
                   "Cell cycle", "TNF signaling pathway",
                   "type I interferon receptor binding", "cellular response to type I interferon",
                   "response to interferon-gamma","locomotion","Jak-STAT signaling pathway", "ncRNA processing")
renamed_terms = c("Antigen presentation", "Cell death", "Cell cycle", "TNF signalling", "Type I IFN", "Type I IFN response", 
                  "IFNg response", "Locomotion", "Jak-STAT signalling", "ncRNA processing")
new_names = data_frame(term.name = selected_terms, new_term_name = renamed_terms)

#Make a plot
plot_dataset = dplyr::filter(go_df, term.name %in% selected_terms) %>% 
  dplyr::mutate(log10p = -log(p.value,10)) %>%
  dplyr::filter(log10p > 8) %>% 
  dplyr::left_join(new_names, by = "term.name")
plot = ggplot(plot_dataset, aes(x = log10p, y = new_term_name)) + 
  facet_grid(cluster_id~.,scales = "free_y", space = "free_y") + 
  geom_point() + 
  scale_x_continuous(limits = c(0,38)) + 
  theme_light() + 
  xlab(expression(paste(-Log[10], " p-value", sep = ""))) +
  theme(strip.text.y = element_text(colour = "grey10"), strip.background = element_rect(fill = "grey85")) +
  theme(axis.title.y = element_blank())
ggsave("figures/main_figures/DE_clusters_GO_enrichment.pdf", plot = plot, width = 2.5, height = 4)



