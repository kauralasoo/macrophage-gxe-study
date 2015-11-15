#Import expression data
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")

#Construct sample metadata for atac
line_metadata = readRDS("macrophage-gxe-study/data/covariates/compiled_line_metadata.rds") #Line metadata
donor_geno_map = dplyr::select(line_metadata, donor, genotype_id) %>% unique()
atac_sample_meta = dplyr::left_join(atac_list$design, donor_geno_map, by = "donor") %>% 
  dplyr::mutate(condition_name = factor(condition, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))
atac_peak_meta = dplyr::mutate(atac_list$gene_metadata, gene_name = gene_id)

#Import genotypes
vcf_file = readRDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")

#### Focus on the CTSC locus ####
#Make a QTL plot for the gene
ctsc_plot = plotEQTL("ENSG00000109861", "rs11019479", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
                     eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)
ggsave("results/SL1344/eQTLs/example_loci/CTSC/CTSC_eQTL.pdf", ctsc_plot, width = 7, height = 7)

#Make a QTL plot for the ATAC peak
ctsc_enhancer = plotEQTL("ATAC_peak_55686", "rs11019479", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_peak_meta)
ggsave("results/SL1344/eQTLs/example_loci/CTSC/CTSC_ATAC_enhancer.pdf", ctsc_enhancer, width = 7, height = 7)

#Downstream and upstream enhancers
ctsc_enhancer1 = plotEQTL("ATAC_peak_55687", "rs11019479", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_peak_meta)
ggsave("results/SL1344/eQTLs/example_loci/CTSC/CTSC_ATAC_enhancer+1.pdf", ctsc_enhancer1, width = 7, height = 7)

ctsc_enhancer_u1 = plotEQTL("ATAC_peak_55685", "rs11019479", atac_list$exprs_cqn, vcf_file$genotypes, 
                          atac_sample_meta, atac_peak_meta)
ggsave("results/SL1344/eQTLs/example_loci/CTSC/CTSC_ATAC_enhancer-1.pdf", ctsc_enhancer_u1, width = 7, height = 7)

#Imprt raw p-values
IFNg_pvalues = readr::read_delim("results/SL1344/fastqtl/output/IFNg_pvalues.txt.gz", delim = " ")

