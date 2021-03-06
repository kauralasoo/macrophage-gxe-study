library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("ggplot2")
library("tidyr")

#Load expression data
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
fastqtl_callset = readRDS("results/SL1344/fastqtl/output/fastqtl_call_set.rds")
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)
vcf_file = readRDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")

#Load GWAS catalog
gwas_catalog = importGwasCatalog("annotations/gwas_catalog_v1.0.1-downloaded_2016-03-02.tsv") %>%
  dplyr::filter(is_european == TRUE, sample_size > 1000)
#Count associations per trait
assocation_count = group_by(gwas_catalog, mapped_trait) %>%
  summarize(association_count = length(mapped_trait))
gwas_catalog_df = dplyr::left_join(gwas_catalog, assocation_count, by = "mapped_trait") %>% 
  dplyr::filter(association_count > 10)

#Import gwas p-values
naive_pvalues = importGwasPvalue("results/SL1344/fastqtl/output/naive_pvalues.coords.txt.gz", gwas_catalog_df)
naive_enrichment = calculateEnrichment(gwas_catalog_df, naive_pvalues, fastqtl_callset$naive, 
                                       gene_id_name_map, p_diff_thresh = 2)

IFNg_pvalues = importGwasPvalue("results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz", gwas_catalog_df)
IFNg_enrichment = calculateEnrichment(gwas_catalog_df, IFNg_pvalues, fastqtl_callset$IFNg,
                                      gene_id_name_map, p_diff_thresh = 2)

SL1344_pvalues = importGwasPvalue("results/SL1344/fastqtl/output/SL1344_pvalues.coords.txt.gz", gwas_catalog_df)
SL1344_enrichment = calculateEnrichment(gwas_catalog_df, SL1344_pvalues, fastqtl_callset$SL1344, 
                                        gene_id_name_map, p_diff_thresh = 2)

IFNg_SL1344_pvalues = importGwasPvalue("results/SL1344/fastqtl/output/IFNg_SL1344_pvalues.coords.txt.gz", gwas_catalog_df)
IFNg_SL1344_enrichment = calculateEnrichment(gwas_catalog_df, IFNg_SL1344_pvalues, fastqtl_callset$IFNg_SL1344, 
                                             gene_id_name_map, p_diff_thresh = 2)

write.table(naive_enrichment$enrichment, "results/SL1344/fastqtl/output/naive_enrichment.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(IFNg_enrichment$enrichment, "results/SL1344/fastqtl/output/IFNg_enrichment.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(SL1344_enrichment$enrichment, "results/SL1344/fastqtl/output/SL1344_enrichment.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(IFNg_SL1344_enrichment$enrichment, "results/SL1344/fastqtl/output/IFNg_SL1344_enrichment.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#Identify all genes that overlap disease loci
overlap_genes = rbind(naive_enrichment$overlap, IFNg_enrichment$overlap, 
                      SL1344_enrichment$overlap, IFNg_SL1344_enrichment$overlap) %>%
  dplyr::arrange(qvalue)
saveRDS(overlap_genes, "results/SL1344/fastqtl/output/disease_overlaps.rds")

#Identify unique set of gene-SNP pairs
disease_overlaps = dplyr::select(overlap_genes, gene_id, gene_name, snp_id, qvalue) %>% 
  dplyr::arrange(qvalue) %>% unique() 
interaction_disease_overlaps = dplyr::semi_join(disease_overlaps, fastqtl_callset$interaction, by = "gene_id") %>% 
  arrange(qvalue)

#Make plots and save them to disk
disease_plots = makeMultiplePlots(interaction_disease_overlaps, eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
                                      eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)
savePlots(disease_plots, "results/SL1344/eQTLs/disease_interaction_plots//", width = 7, height = 7)

#Combine all enrichments
dplyr::select(naive_pvalues, mapped_trait, OR_trait_qtls)

naive_OR = dplyr::transmute(naive_enrichment$enrichment,mapped_trait, naive = OR_trait_qtls)
ifng_OR = dplyr::transmute(IFNg_enrichment$enrichment,mapped_trait, IFNg = OR_trait_qtls)
sl1344_OR = dplyr::transmute(SL1344_enrichment$enrichment,mapped_trait, SL1344 = OR_trait_qtls)
ifng_sl1344_OR = dplyr::transmute(IFNg_SL1344_enrichment$enrichment,mapped_trait, IFNg_SL1344 = OR_trait_qtls)

enrichment_table = dplyr::left_join(naive_OR, ifng_OR, by = "mapped_trait") %>% 
  dplyr::left_join(sl1344_OR, by = "mapped_trait") %>% 
  dplyr::left_join(ifng_sl1344_OR, by ="mapped_trait") %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(max_OR = max(naive, IFNg, SL1344, IFNg_SL1344)) %>% 
  ungroup() %>% 
  dplyr::arrange(-max_OR) %>%
  dplyr::top_n(10) %>%
  dplyr::mutate(mapped_trait = factor(mapped_trait, rev(mapped_trait))) %>%
  tidyr::gather(key = condition, value = OR, naive:IFNg_SL1344)
  
trait_enrichment_plot = ggplot(enrichment_table, aes(x = OR, y = mapped_trait, color = condition)) + 
  geom_point()
ggsave("results/SL1344/eQTLs/disease_enrichments.pdf", trait_enrichment_plot, width = 8, height = 5)

