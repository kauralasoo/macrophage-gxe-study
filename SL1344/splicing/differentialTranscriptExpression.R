library("DRIMSeq")
library("SummarizedExperiment")

se_ensembl = readRDS("results/SL1344/combined_ensembl85_transcript_quants.rds")

#Construct splicing dataset
splicing_data <- DRIMSeq::dmDSdata(counts = assays(se_ensembl)$counts, 
              gene_id = rowData(se_ensembl)$gene_id,
              feature_id = rowData(se_ensembl)$transcript_id, 
              sample_id = se_ensembl$sample_id, 
              group = se_ensembl$condition_name)

#Filter by counts
splicing_data <- dmFilter(splicing_data, min_samps_gene_expr = 84, min_samps_feature_expr = 84,
              min_samps_feature_prop = 0)

#Select a subset of genes
#subset_genes = c("ENSG00000111912", "ENSG00000151715", "ENSG00000184517")
#splicing_data <- splicing_data[names(splicing_data) %in% subset_genes, ]

#Estimate dispersions
splicing_data <- dmDispersion(splicing_data, verbose = 1, 
                              BPPARAM = BiocParallel::SerialParam())
#Save dispersion estimates to disk
saveRDS(splicing_data, "results/SL1344/combined_ensembl85_dispersion_estimates.serial.rds")
splicing_data = readRDS("results/SL1344/combined_ensembl85_dispersion_estimates.serial.rds")

#Estimate proportions
splicing_data <- dmFit(splicing_data, BPPARAM = BiocParallel::SerialParam())

#Test for differences in proportions
all_dte <- dmTest(splicing_data, verbose = 1, BPPARAM = BiocParallel::SerialParam())
naive_vs_SL1344_dte <- dmTest(splicing_data, verbose = 1, BPPARAM = BiocParallel::SerialParam(),
                        compared_groups = c(1,3))
naive_vs_ifng_dte <- dmTest(splicing_data, verbose = 1, BPPARAM = BiocParallel::SerialParam(),
                              compared_groups = c(1,2))
naive_vs_ifng_sl1344_dte <- dmTest(splicing_data, verbose = 1, BPPARAM = BiocParallel::SerialParam(),
                            compared_groups = c(1,4))

