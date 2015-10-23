#Compress and index SNP coordinates file
bgzip results/SL1344/eqtlbma/input/snp_coords.bed 
tabix -p bed results/SL1344/eqtlbma/input/snp_coords.bed.gz

#Split genes into batches
/software/R-3.1.2/bin/Rscript ~/software/utils/eqtlbma/splitGeneCoords.R -i results/SL1344/eqtlbma/input/gene_coords.bed -b 200 -o results/SL1344/eqtlbma/input/gene_batches.txt

#Compute bayes factors for each SNP/gene pair
cat results/SL1344/eqtlbma/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 6000 --jobname run_eqtlbma_bf --command  "python ~/software/utils/eqtlbma/run_eqtlbma_bf.py --indir results/SL1344/eqtlbma/input/ --outdir results/SL1344/eqtlbma/output/ --outprefix eqtlbma"

#Run eqtlbma with permutations
#Run each batch with 250 permutations
cat results/SL1344/eqtlbma/input/gene_batches.txt | tail -n 82 | python ~/software/utils/submitJobs.py --MEM 15000 --jobname run_eqtlbma_bf_perm --ncores 4 --command  "python ~/software/utils/eqtlbma/run_eqtlbma_bf.py --indir results/SL1344/eqtlbma/input/ --outdir results/SL1344/eqtlbma/output_perm/ --outprefix eqtlbma_perm --thread 4 --nperm 250"

#MHC genes take longer to run
grep chr_6_batch_1 results/SL1344/eqtlbma/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 44000 --jobname run_eqtlbma_bf_perm --ncores 4 --queue hugemem --command  "python ~/software/utils/eqtlbma/run_eqtlbma_bf.py --indir results/SL1344/eqtlbma/input/ --outdir results/SL1344/eqtlbma/output_perm/ --outprefix eqtlbma_perm --thread 4 --nperm 250"
grep chr_6_batch_2 results/SL1344/eqtlbma/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 44000 --jobname run_eqtlbma_bf_perm --ncores 4 --queue hugemem --command  "python ~/software/utils/eqtlbma/run_eqtlbma_bf.py --indir results/SL1344/eqtlbma/input/ --outdir results/SL1344/eqtlbma/output_perm/ --outprefix eqtlbma_perm --thread 4 --nperm 250"

#Merge results from all batches
python ~/software/utils/eqtlbma/merge_eqtlbma_bf.py --geneBatches results/SL1344/eqtlbma/input/gene_batches.txt --outprefix eqtlbma --outdir results/SL1344/eqtlbma/output/

#Run the hierarchical model
echo "None" | python ~/software/utils/submitJobs.py --MEM 15000 --jobname run_eqtlbma_hm --ncores 6 --command "python ~/software/utils/eqtlbma/run_eqtlbma_hm.py --workdir results/SL1344/eqtlbma/output/ --outprefix eqtlbma --nsubgrp 4 --dim 15 --thread 6"

#Run the EBF procedure
bsub -G team170 -n1 -R "span[hosts=1] select[mem>12000] rusage[mem=12000]" -q normal -M 12000 -o FarmOut/eqtlbma_EBF.%J.jobout "python ~/software/utils/eqtlbma/run_eqtlbma_avg_bfs.py --workdir eqtlbma/output_imputed/ --outprefix eqtlbma_imputed --nsubgrp 4 --dim 15 --mode ebf"

#Calculate average BFS over configurations
bsub -G team170 -n1 -R "span[hosts=1] select[mem>12000] rusage[mem=12000]" -q normal -M 12000 -o FarmOut/eqtlbma_avg_bfs.%J.jobout "python ~/software/utils/eqtlbma/run_eqtlbma_avg_bfs.py --workdir eqtlbma/output_imputed/ --outprefix eqtlbma_imputed --nsubgrp 4 --dim 15 --mode post --pi0 0.6165"
