library("devtools")
library("qvalue")
library("dplyr")
library("ggplot2")
library("SNPRelate")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")
library("DESeq2")

#Import expression data
eqtl_data_list = readRDS("results/acLDL/acLDL_eqtl_data_list.rds")
expression_data_list = readRDS("results/acLDL/acLDL_combined_expression_data.rds")
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Import genotypes
vcf_file = readRDS("results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.rds")

#Extract counts
total_counts = expression_data_list$exprs_counts[,eqtl_data_list$sample_metadata$sample_id]

#Extract LIPA genotype
LIPA_genotype = vcf_file$genotypes["rs1412444",]
LIPA_geno_df = data_frame(genotype_id = names(LIPA_genotype), rs1412444 = LIPA_genotype)
design = dplyr::left_join(eqtl_data_list$sample_metadata, LIPA_geno_df, by = "genotype_id") %>% as.data.frame()
rownames(design) = design$sample_id

dds = DESeqDataSetFromMatrix(total_counts, design, ~condition + rs1412444 + condition:rs1412444)
dds = DESeq(dds)

stimulation_res = results(dds, name = "condition_Ctrl_vs_AcLDL")
stimulation_results = filterDESeqResults(stimulation_res,eqtl_data_list$gene_metadata, min_padj = 0.01, min_fc = 0)

snp_res = results(dds, name = "rs1412444")
snp_results = filterDESeqResults(snp_res,eqtl_data_list$gene_metadata, min_padj = 0.01, min_fc = 0)

inter_res = results(dds, name = "conditionCtrl.rs1412444")
inter_results = filterDESeqResults(snp_res,eqtl_data_list$gene_metadata, min_padj = 0.01, min_fc = 0)

