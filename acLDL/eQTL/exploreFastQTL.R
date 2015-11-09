library("devtools")
library("qvalue")
library("dplyr")
library("ggplot2")
library("SNPRelate")
load_all("../seqUtils/")
load_all("macrophage-gxe-study/housekeeping/")

#Import expression data
eqtl_data_list = readRDS("results/acLDL/acLDL_eqtl_data_list.rds")
gene_id_name_map = dplyr::select(eqtl_data_list$gene_metadata, gene_id, gene_name)

#Import genotypes
SNPRelate::snpgdsVCF2GDS("results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz", 
"results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.gds", method = "copy.num.of.ref")
vcf_file = gdsToMatrix("results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.gds")
saveRDS(vcf_file, "results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.rds")
vcf_file = readRDS("results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.rds")


#Import permutation eQTL results
ctrl_qtls = importFastQTLTable("results/acLDL/fastqtl/output/Ctrl_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)
acLDL_qtls = importFastQTLTable("results/acLDL/fastqtl/output/AcLDL_permuted.txt.gz") %>% enrichFastQTLPvalues(gene_id_name_map)

#Calculate Pi1 between stimulated and unstimulated
qtl_list = list(Ctrl = ctrl_qtls, AcLDL = acLDL_qtls)
pi1_matrix = calculatePairwisePi1(qtl_list)

#Replication between the two studies is relatively low, BUT they individuals and variants tested do not overlap completely.

#Find hits
ctrl_hits = dplyr::filter(ctrl_qtls, qvalue < 0.1)
acldl_hits = dplyr::filter(acLDL_qtls, qvalue < 0.1)
joint_qtls = rbind(ctrl_hits, acldl_hits) %>% dplyr::select(gene_id, snp_id) %>% unique()

interactions = testMultipleInteractions(joint_qtls, eqtl_data_list, vcf_file)
interaction_df = plyr::ldply(interactions, .id = "id")
saveRDS(interactions, "results/acLDL/fastqtl/output/interactions.rds")

#Make a plot of the LIPA qtl
lipa_qtl = plotEQTL("ENSG00000107798", "rs1332329", eqtl_data_list$exprs_cqn, vcf_file$genotypes,
         eqtl_data_list$sample_metadata, eqtl_data_list$gene_metadata)
ggsave("results/acLDL/LIPA_qtl.pdf", lipa_qtl, width = 8, height = 6)
