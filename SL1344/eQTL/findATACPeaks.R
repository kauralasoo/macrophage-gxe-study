library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("MatrixEQTL")

#Import expression data
eqtl_data_list = readRDS("results/SL1344/eqtl_data_list.rds")
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Import ATAC data
atac_list = readRDS("../macrophage-chromatin/results/ATAC/ATAC_combined_accessibility_data.rds")

#Import genotypes
vcf_file = readRDS("results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.rds")

#Import eQTL mapping hits
fastqtl_callset = readRDS("results/SL1344/fastqtl/output/fastqtl_call_set.rds")
interaction_snps = unique(fastqtl_callset$interaction$snp_id)

#Construct genepos and SNP-pos matrices
atac_snpspos = dplyr::filter(vcf_file$snpspos, snpid %in% interaction_snps) %>% 
  as.data.frame()
atac_genepos = dplyr::rename(atac_list$gene_metadata, geneid = gene_id) %>% 
  dplyr::select(geneid, chr, left, right)
IFNg_atac_hits = runMatrixEQTL(exp_data = atac_list$exprs_cqn_list$IFNg, 
                    geno_data = extractSubset(atac_list$design_list$IFNg, vcf_file$genotypes, "genotype_id", "donor"),
                    snpspos = atac_snpspos, 
                    genepos = atac_genepos, 
                    cisDist = 2e05, pvOutputThreshold = 1)

naive_atac_hits = runMatrixEQTL(exp_data = atac_list$exprs_cqn_list$naive, 
                               geno_data = extractSubset(atac_list$design_list$naive, vcf_file$genotypes, "genotype_id", "donor"),
                               snpspos = atac_snpspos, 
                               genepos = atac_genepos, 
                               cisDist = 2e05, pvOutputThreshold = 1)



