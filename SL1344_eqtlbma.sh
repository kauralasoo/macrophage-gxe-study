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