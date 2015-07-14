labels = c("gene_id","snp_id","chromosome", "snp_pos", "REF","ALT","AAF","HWE","IA",
           "R2","chi_sq","Pi","Phi","Delta", "OD","snp_id_region","fSNP_count",
           "tSNP_count","null_iter","alt_iter","ties","null_lik","convergence",
           "fSNP_corr","rSNP_corr")
cd14 = read.table("rasqual/input/CD14.txt")
colnames(cd14) = labels
plot(cd14$snp_pos,cd14$chi_sq)