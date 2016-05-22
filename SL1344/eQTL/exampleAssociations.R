load_all("../seqUtils/")

#Import expression data
combined_expression_data = readRDS("results/SL1344/combined_expression_data_covariates.rds")
gene_id_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)
combined_expression_data$sample_metadata$condition_name = factor(combined_expression_data$sample_metadata$condition_name, 
        levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data_covariates.rds")
atac_list$sample_metadata = dplyr::mutate(atac_list$sample_metadata, 
          condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344")))

#Import genotypes
vcf_file = readRDS("genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.rds")


#### Focus on the CTSC locus ####
#Make a QTL plot for the gene
ctsc_plot = plotEQTL("ENSG00000144228", "rs12621644", combined_expression_data$cqn, vcf_file$genotypes, 
                     combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)
ggsave("results/SL1344/eQTLs/example_loci/CTSC/CTSC_eQTL.pdf", ctsc_plot, width = 7, height = 7)

#Make a QTL plot for the ATAC peak
ctsc_enhancer = plotEQTL("ATAC_peak_55686", "rs11019479", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata)
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

#GP1BA
gp1ba_enhancer = plotEQTL("ATAC_peak_113648", "rs238242", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata)
ggsave("results/SL1344/eQTLs/example_loci/GP1BA/GP1BA_ATAC_enhancer.pdf", gp1ba_enhancer, width = 7, height = 7)

plotEQTL("ENSG00000185245", "rs238242", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
                     eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata) %>%
ggsave("results/SL1344/eQTLs/example_loci/GP1BA/GP1BA_eQTL.pdf", ., width = 7, height = 7)


#SPOPL
plotEQTL("ATAC_peak_154955", "rs12621644", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/SPOPL/SPOPL_enhancer.pdf", ., width = 7, height = 7)

plotEQTL("ATAC_peak_154957", "rs12621644", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/SPOPL/SPOPL_enhancer+2.pdf", ., width = 7, height = 7)

plotEQTL("ATAC_peak_154958", "rs12621644", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/SPOPL/SPOPL_enhancer+3.pdf", ., width = 7, height = 7)

spopl_data = plotEQTL("ENSG00000144228", "rs12621644", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata, return_df = TRUE) 
plotQtlRow(spopl_data) %>% 
  ggsave("results/SL1344/eQTLs/example_loci/SPOPL/SPOPL_eQTL.pdf", ., width = 6, height = 2.5)




#CCS
plotEQTL("ATAC_peak_52919", "rs566673", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/CCS/CCS_enhancer.pdf", ., width = 7, height = 7)

plotEQTL("ENSG00000173992", "rs566673", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
         eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/CCS/CCS_eQTL.pdf", ., width = 7, height = 7)

#RGS14
plotEQTL("ATAC_peak_239461", "rs12654812", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/RGS14/RGS14_enhancer1.pdf", ., width = 7, height = 7)

plotEQTL("ATAC_peak_239459", "rs12654812", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/RGS14/RGS14_enhancer2.pdf", ., width = 7, height = 7)

plotEQTL("ATAC_peak_239457", "rs12654812", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/RGS14/RGS14_enhancer2-2.pdf", ., width = 7, height = 7)

plotEQTL("ENSG00000169220", "rs12654812",combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata, return_df = TRUE) %>%
  plotQtlRow() %>% 
  ggsave("results/SL1344/eQTLs/example_loci/RGS14/RGS14_eQTL.pdf", ., width = 6, height = 2.5)

plotEQTL("ENSG00000146094", "rs12654812", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
         eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)  %>%
  ggsave("results/SL1344/eQTLs/example_loci/RGS14/DOK3_eQTL.pdf", ., width = 7, height = 7)

#TMEM229B
plotEQTL("ATAC_peak_88584", "rs12878807", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/TMEM229B/TMEM229B_enhancer.pdf", ., width = 7, height = 7)

plotEQTL("ATAC_peak_88585", "rs12878807", atac_list$exprs_cqn, vcf_file$genotypes, 
         atac_sample_meta, atac_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/TMEM229B/TMEM229B_enhancer+1.pdf", ., width = 7, height = 7)


plotEQTL("ENSG00000198133", "rs12878807", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
         eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/TMEM229B/TMEM229B_eQTL.pdf", ., width = 7, height = 7)

plotEQTL("ENSG00000054690", "rs12878807", eqtl_data_list$exprs_cqn, vcf_file$genotypes, 
         eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata) %>%
  ggsave("results/SL1344/eQTLs/example_loci/TMEM229B/PLEKHH1_eQTL.pdf", ., width = 7, height = 7)


plotEQTL("ENSG00000197728", "rs1131017", combined_expression_data$cqn, vcf_file$genotypes, 
         combined_expression_data$sample_metadata, combined_expression_data$gene_metadata)