library("dplyr")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")

full_pvalues = list(CD14 = importQTLtoolsTable("processed/Fairfax/qtltools/output/full/CD14.permuted.txt.gz"),
                    IFN = importQTLtoolsTable("processed/Fairfax/qtltools/output/full/IFN.permuted.txt.gz"),
                    LPS2 = importQTLtoolsTable("processed/Fairfax/qtltools/output/full/LPS2.permuted.txt.gz"),
                    LPS24 = importQTLtoolsTable("processed/Fairfax/qtltools/output/full/LPS24.permuted.txt.gz"))
shared_pvalues = list(CD14 = importQTLtoolsTable("processed/Fairfax/qtltools/output/shared/CD14.permuted.txt.gz"),
                    IFN = importQTLtoolsTable("processed/Fairfax/qtltools/output/shared/IFN.permuted.txt.gz"),
                    LPS2 = importQTLtoolsTable("processed/Fairfax/qtltools/output/shared/LPS2.permuted.txt.gz"),
                    LPS24 = importQTLtoolsTable("processed/Fairfax/qtltools/output/shared/LPS24.permuted.txt.gz"))
shared_84_pvalues = list(CD14 = importQTLtoolsTable("processed/Fairfax/qtltools/output/shared_84/CD14.permuted.txt.gz"),
                    IFN = importQTLtoolsTable("processed/Fairfax/qtltools/output/shared_84/IFN.permuted.txt.gz"),
                    LPS2 = importQTLtoolsTable("processed/Fairfax/qtltools/output/shared_84/LPS2.permuted.txt.gz"),
                    LPS24 = importQTLtoolsTable("processed/Fairfax/qtltools/output/shared_84/LPS24.permuted.txt.gz"))

#Put all results into a list
QTL_min_pvalue_list = list(full = full_pvalues, shared = shared_pvalues, shared_84 = shared_84_pvalues)
saveRDS(QTL_min_pvalue_list, "results/Fairfax/fairfax_qtl_min_pvalues.rds")
