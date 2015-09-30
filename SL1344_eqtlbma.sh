#Remove SNPs with multiple IDs
grep -v ";" eqtlbma/snp_coords.bed > eqtlbma/snp_coords_filtered.bed
cat eqtlbma/snp_coords_filtered.bed | sort -k1,1V -k2,2g | bgzip > eqtlbma/snp_coords_filtered.bed.gz
tabix -p bed eqtlbma/snp_coords_filtered.bed.gz 

#Run eqtlbma on a small number of genes
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o FarmOut/eqtlbma_bf.%J.jobout "eqtlbma_bf --geno eqtlbma/list_genotypes.txt --scoord eqtlbma/snp_coords.bed.gz --exp eqtlbma/list_explevels.txt --gcoord eqtlbma/gene_coords.subset.sort.bed --anchor TSS --cis 500000 --out eqtlbma/output/subset_test --analys join --gridL eqtlbma/grid_phi2_oma2_general.txt  --gridS eqtlbma/grid_phi2_oma2_with-configs.txt --bfs all --error mvlr -v 3"

bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o FarmOut/eqtlbma_hm.%J.jobout  "eqtlbma_hm --data 'eqtlbma/output/subset_test_l10abfs_raw.txt.gz' --nsubgrp 4 --dim 15 --ngrid 10 --out eqtlbma/output/subset_test_eqtlbma_hm_txt.gz"

#Extract config weights
zcat eqtlbma/output/subset_test_eqtlbma_hm_txt.gz | grep "#config" | awk '{split($1,a,"."); print a[2]"\t"$2}' > eqtlbma/output/subset_test_config_weights.txt

#Extract grid weigths
zcat eqtlbma/output/subset_test_eqtlbma_hm_txt.gz | grep "#grid" | cut -f2 > eqtlbma/output/subset_test_grid_weights.txt

eqtlbma_avg_bfs --in eqtlbma/output/subset_test_l10abfs_raw.txt.gz --gwts eqtlbma/output/subset_test_grid_weights.txt --nsubgrp 4 --dim 15 --cwts eqtlbma/output/subset_test_config_weights.txt --save bf+post --pi0 0.783629 --post a+b+c+d --bestdim --alldim --out eqtlbma/output/out_eqtlbma_avg_bfs.txt.gz --thread 1


#Run eqtlbma with covariates
#Compute bayes factors for each SNP/gene pair
bsub -G team170 -n1 -R "span[hosts=1] select[mem>4000] rusage[mem=4000]" -q normal -M 4000 -o FarmOut/eqtlbma_bf.%J.jobout "eqtlbma_bf --geno eqtlbma/input/list.genotypes.txt --scoord eqtlbma/input/snp_coords.bed.gz --exp eqtlbma/input/list.expression.txt --gcoord eqtlbma/input/gene_coords.bed --anchor TSS --cis 500000 --out eqtlbma/output/out_eqtlbma_full --analys join --covar eqtlbma/input/list.covariates.txt --gridL macrophage-gxe-study/data/eqtlbma/grid_phi2_oma2_general.txt  --gridS macrophage-gxe-study/data/eqtlbma/grid_phi2_oma2_with-configs.txt --bfs all --error mvlr -v 3"

#Run the hierarchical model
bsub -G team170 -n1 -R "span[hosts=1] select[mem>6000] rusage[mem=6000]" -q normal -M 6000 -o FarmOut/eqtlbma_hm.%J.jobout "python ~/software/utils/run_eqtlbma_hm.py --workdir eqtlbma/output/ --outprefix out_eqtlbma_full --nsubgrp 4 --dim 15"

#Run the EBF procedure
bsub -G team170 -n1 -R "span[hosts=1] select[mem>8000] rusage[mem=8000]" -q normal -M 8000 -o FarmOut/eqtlbma_EBF.%J.jobout "eqtlbma_avg_bfs --in 'eqtlbma/output/out_eqtlbma_full_l10abfs_raw.txt.gz' --gwts eqtlbma/output/out_eqtlbma_full_grid_weights.txt --nsubgrp 4 --dim 15 --cwts eqtlbma/output/out_eqtlbma_full_config_weights.txt --save bf --out eqtlbma/output/out_eqtlbma_avg_bfs_genes.txt.gz"

#Calculate average BFS over configurations
bsub -G team170 -n1 -R "span[hosts=1] select[mem>8000] rusage[mem=8000]" -q normal -M 8000 -o FarmOut/eqtlbma_avg_bfs.%J.jobout "eqtlbma_avg_bfs --in 'eqtlbma/output/out_eqtlbma_full_l10abfs_raw.txt.gz' --gwts eqtlbma/output/out_eqtlbma_full_grid_weights.txt --nsubgrp 4 --dim 15 --cwts eqtlbma/output/out_eqtlbma_full_config_weights.txt --save bf+post --pi0 0.6673715 --post a+b+c+d --bestdim --alldim --out eqtlbma/output/out_eqtlbma_avg_bfs.txt.gz --thread 1"