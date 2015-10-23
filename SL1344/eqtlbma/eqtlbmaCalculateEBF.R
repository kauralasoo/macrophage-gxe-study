source("~/farm-home/software/eqtlbma/scripts/utils_eqtlbma.R")

#Run the EBF procedure
gene.bfs <- read.table("results/SL1344/eqtlbma/output/eqtlbma_avg_bfs_EBF.txt.gz", header=TRUE)
pi0.ebf <- estimatePi0WithEbf(log10.bfs=gene.bfs$gene.log10.bf[!duplicated(gene.bfs$gene)],verbose=1)
called.nulls.ebf <- ! controlBayesFdr(log10.bfs=gene.bfs$gene.log10.bf[!duplicated(gene.bfs$gene)],
                                      pi0=pi0.ebf, fdr.level=0.1, verbose=1)

#Run the QBF procedure
gene.bfs <- read.table("results/SL1344/eqtlbma/output/eqtlbma_avg_bfs_EBF.txt.gz", header=TRUE)
perms.qbf = read.table("results/SL1344/eqtlbma/output_perm/eqtlbma_perm_joinPermPvals.txt.gz", header = TRUE)
pi0.qbf <- estimatePi0WithQbf(log10.bfs=as.matrix(perms.qbf[!duplicated(perms.qbf$gene),c(5,6)]),gamma=0.5, verbose=1)
called.nulls.qbf <- ! controlBayesFdr(log10.bfs=gene.bfs$gene.log10.bf[!duplicated(gene.bfs$gene)],
                                      pi0=pi0.qbf, fdr.level=0.1, verbose=1)