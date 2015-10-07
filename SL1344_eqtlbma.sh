#Compress and index SNP coordinates file
bgzip eqtlbma/input_imputed/snp_coords.bed 
tabix -p bed eqtlbma/input_imputed/snp_coords.bed.gz

#Split genes into batches
/software/R-3.1.2/bin/Rscript ~/software/utils/eqtlbma/splitGeneCoords.R -i eqtlbma/input_imputed/gene_coords.bed -b 200 -o eqtlbma/input_imputed/gene_batches.txt

#Compute bayes factors for each SNP/gene pair
cat eqtlbma/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname run_eqtlbma_bf --command  "python ~/software/utils/eqtlbma/run_eqtlbma_bf.py --indir eqtlbma/input_imputed/ --outdir eqtlbma/output_imputed/ --outprefix eqtlbma_imputed"

#Merge results from all batches
python ~/software/utils/eqtlbma/merge_eqtlbma_bf.py --geneBatches eqtlbma/input_imputed/gene_batches.txt --outprefix eqtlbma_imputed --outdir eqtlbma/output_imputed/

#Run the hierarchical model
echo "None" | python ~/software/utils/submitJobs.py --MEM 12000 --jobname run_eqtlbma_hm --ncores 6 --command "python ~/software/utils/eqtlbma/run_eqtlbma_hm.py --workdir eqtlbma/output_imputed/ --outprefix eqtlbma_imputed --nsubgrp 4 --dim 15 --thread 6"

#Keep only a subset of the configurations
bsub -G team170 -n1 -R "span[hosts=1] select[mem>12000] rusage[mem=12000]" -q normal -M 12000 -o FarmOut/eqtlbma_hm.%J.jobout "python ~/software/utils/eqtlbma/run_eqtlbma_hm.py --workdir eqtlbma/output_/ --outprefix eqtlbma_configs --nsubgrp 4 --dim 15 --configs '1-2-3-4|1-2-3|1-2|2-3|1-4|2-4'"

#Run the EBF procedure
bsub -G team170 -n1 -R "span[hosts=1] select[mem>12000] rusage[mem=12000]" -q normal -M 12000 -o FarmOut/eqtlbma_EBF.%J.jobout "python ~/software/utils/eqtlbma/run_eqtlbma_avg_bfs.py --workdir eqtlbma/output_imputed/ --outprefix eqtlbma_imputed --nsubgrp 4 --dim 15 --mode ebf"

#Calculate average BFS over configurations
bsub -G team170 -n1 -R "span[hosts=1] select[mem>12000] rusage[mem=12000]" -q normal -M 12000 -o FarmOut/eqtlbma_avg_bfs.%J.jobout "python ~/software/utils/eqtlbma/run_eqtlbma_avg_bfs.py --workdir eqtlbma/output_imputed/ --outprefix eqtlbma_imputed --nsubgrp 4 --dim 15 --mode post --pi0 0.6417"
