source("~/farm-home/software/eqtlbma/scripts/utils_eqtlbma.R")

gene.bfs <- read.table("eqtlbma/output_imputed/eqtlbma_imputed_avg_bfs_EBF.txt.gz", header=TRUE)
pi0.ebf <- estimatePi0WithEbf(log10.bfs=gene.bfs$gene.log10.bf[!duplicated(gene.bfs$gene)],verbose=1)
called.nulls.ebf <- ! controlBayesFdr(log10.bfs=gene.bfs$gene.log10.bf[!duplicated(gene.bfs$gene)],
                                      pi0=pi0.ebf, fdr.level=0.1, verbose=1)