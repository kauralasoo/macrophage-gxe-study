load_all("../seqUtils/")
library("ggplot2")

#Import p-values
ensembl_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/AcLDL.permuted.txt.gz"),
  Diff = importQTLtoolsTable("processed/acLDL/fastqtl_output/ensembl_87/Diff.permuted.txt.gz"))
revised_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/reviseAnnotations/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/fastqtl_output/reviseAnnotations/AcLDL.permuted.txt.gz"),
  Diff = importQTLtoolsTable("processed/acLDL/fastqtl_output/reviseAnnotations/Diff.permuted.txt.gz"))
leafcutter_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/fastqtl_output/leafcutter/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/fastqtl_output/leafcutter/AcLDL.permuted.txt.gz"),
  Diff = importQTLtoolsTable("processed/acLDL/fastqtl_output/leafcutter/Diff.permuted.txt.gz"))


#Make qq-plots
qq_df = dplyr::mutate(ensembl_pvalues$Ctrl, p_eigen = p_perm) %>% dplyr::arrange(p_eigen) %>% addExpectedPvalue()
ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_eigen,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")