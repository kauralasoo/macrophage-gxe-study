library("dplyr")
library("devtools")
load_all("../seqUtils/")

#Load expression data
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
fastqtl_callset = readRDS("results/SL1344/fastqtl/output/fastqtl_call_set.rds")

#Load GWAS catalog
gwas_catalog = importGwasCatalog("annotations/gwas_catalog_v1.0.1-downloaded_2015-11-19.tsv") %>%
  dplyr::filter(is_european == TRUE, sample_size > 1000)
#Count associations per trait
assocation_count = group_by(gwas_catalog, mapped_trait) %>%
  summarize(association_count = length(mapped_trait))
gwas_catalog_df = dplyr::left_join(gwas_catalog, assocation_count, by = "mapped_trait") %>% 
  dplyr::filter(association_count > 10)

#Import all p-values
n_pval = readr::read_delim("results/SL1344/fastqtl/output/SL1344_pvalues.coords.txt.gz", delim = " ", col_names = FALSE)
colnames(n_pval) = c("gene_id", "chr","pos", "snp_id", "distance", "pvalue", "effect_size")

#For each trait-SNP pair find the most associated gene
gwas_pvalues = dplyr::semi_join(n_pval, gwas_catalog_df, by = c("chr","pos"))
genes_with_traits = dplyr::inner_join(gwas_pvalues, gwas_catalog_df, by = c("chr", "pos")) %>%
  dplyr::select(gene_id, chr, pos, pvalue, trait, mapped_trait) %>% 
  dplyr::group_by(chr, pos, mapped_trait) %>% 
  dplyr::arrange(pvalue) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(gene_id, chr, pos, pvalue, mapped_trait) %>% 
  dplyr::group_by(gene_id, mapped_trait) %>% 
  dplyr::arrange(pvalue) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::rename(trait_chr = chr, trait_pos = pos, trait_pvalue = pvalue)

#Count the number of idendent genes per trait
genes_per_trait = genes_with_traits %>% 
  group_by(mapped_trait) %>% 
  summarize(trait_gene_count = length(gene_id)) %>% 
  arrange(-trait_gene_count)
  
#Join QTLs and GWAS hits
trait_qtl_map = dplyr::inner_join(fastqtl_callset$SL1344, genes_with_traits, by = "gene_id") %>%
  dplyr::mutate(pvalue_diff = -log10(p_nominal) +log10(trait_pvalue)) %>% 
  dplyr::filter(pvalue_diff < 3) #Difference in eQTL pvalues at most 3 orders of magnitude

#Count the number QTL genes that overlap a trait
trait_qtl_overlap = trait_qtl_map %>% 
  group_by(mapped_trait) %>% 
  dplyr::summarise(trait_qtl_overlap = length(gene_id)) %>% 
  arrange(-trait_qtl_overlap)

#Calculate odd-ratios of enrichment
n_qtls = nrow(fastqtl_callset$naive)
n_genes = nrow(eqtl_data_list$exprs_cqn)
overlap_df = dplyr::left_join(genes_per_trait, trait_qtl_overlap, by = "mapped_trait") %>% 
  dplyr::mutate(trait_qtl_overlap = ifelse(is.na(trait_qtl_overlap), 0, trait_qtl_overlap)) %>%
  dplyr::mutate(OR = (trait_qtl_overlap/(n_qtls - trait_qtl_overlap)) / ((trait_gene_count-trait_qtl_overlap)/(n_genes-n_qtls))) %>%
  dplyr::arrange(-OR) %>%
  dplyr::filter(trait_gene_count > 20)
