#RUN rasqual on all peaks
#naive
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.gc_library_size.bin --n 42 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/naive_100kb/batches/naive_100kb --covariates results/ATAC/rasqual/input/naive.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5 --fix-genotype'"
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/IFNg.expression.bin --offsets results/ATAC/rasqual/input/IFNg.gc_library_size.bin --n 41 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/IFNg_100kb/batches/IFNg_100kb --covariates results/ATAC/rasqual/input/IFNg.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5 --fix-genotype'"
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/SL1344.expression.bin --offsets results/ATAC/rasqual/input/SL1344.gc_library_size.bin --n 31 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/SL1344.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/SL1344_100kb/batches/SL1344_100kb --covariates results/ATAC/rasqual/input/SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5 --fix-genotype'"
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/IFNg_SL1344.expression.bin --offsets results/ATAC/rasqual/input/IFNg_SL1344.gc_library_size.bin --n 31 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches/IFNg_SL1344_100kb --covariates results/ATAC/rasqual/input/IFNg_SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5 --fix-genotype'"

#Merge batches into single files
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/ATAC/rasqual/output/naive_100kb/batches/naive_100kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/ATAC/rasqual/output/IFNg_100kb/batches/IFNg_100kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/ATAC/rasqual/output/SL1344_100kb/batches/SL1344_100kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches/IFNg_SL1344_100kb"

#Move merged files to parent directory
mv results/ATAC/rasqual/output/naive_100kb/batches/naive_100kb.txt results/ATAC/rasqual/output/naive_100kb/
mv results/ATAC/rasqual/output/IFNg_100kb/batches/IFNg_100kb.txt results/ATAC/rasqual/output/IFNg_100kb/
mv results/ATAC/rasqual/output/SL1344_100kb/batches/SL1344_100kb.txt results/ATAC/rasqual/output/SL1344_100kb/
mv results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches/IFNg_SL1344_100kb.txt results/ATAC/rasqual/output/IFNg_SL1344_100kb/

# Compress batches
tar czf results/ATAC/rasqual/output/naive_100kb/batches.tar.gz results/ATAC/rasqual/output/naive_100kb/batches/ &
tar czf results/ATAC/rasqual/output/IFNg_100kb/batches.tar.gz results/ATAC/rasqual/output/IFNg_100kb/batches/ &
tar czf results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches.tar.gz results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches/ &
tar czf results/ATAC/rasqual/output/SL1344_100kb/batches.tar.gz results/ATAC/rasqual/output/SL1344_100kb/batches/ &

#Extract completed gene ids
cut -f1 results/ATAC/rasqual/output/naive_100kb/naive_100kb.txt | uniq > results/ATAC/rasqual/output/naive_100kb/naive_100kb.completed_ids.txt &
cut -f1 results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.txt | uniq > results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.completed_ids.txt &
cut -f1 results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.txt | uniq > results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.completed_ids.txt &
cut -f1 results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.txt | uniq > results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.completed_ids.txt &

#Sort and index rasqual output for tabix
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED  results/ATAC/rasqual/output/naive_100kb/naive_100kb.txt | sort -k3,3 -k4,4n | bgzip > results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.txt | sort -k3,3 -k4,4n | bgzip > results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout  "grep -v SKIPPED  results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.txt | sort -k3,3 -k4,4n | bgzip > results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout  "grep -v SKIPPED  results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.txt | sort -k3,3 -k4,4n | bgzip > results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz"

#Convert rasqual output into format suitable for eigenMT
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/ATAC/rasqual/output/naive_100kb/naive_100kb.txt > results/ATAC/rasqual/output/naive_100kb/naive_100kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.txt > results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.txt > results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.txt > results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.eigenMT_input.txt"

#Run eigenMT chromosome-by-chromosme
#Make sure that peak locations are properly exported
cat ../macrophage-gxe-study/macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 4 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir ../macrophage-gxe-study/results/SL1344/eigenMT/input/ --genepos results/ATAC/eigenMT/input/gene_positions.txt --QTL results/ATAC/rasqual/output/naive_100kb/naive_100kb.eigenMT_input.txt --out_prefix results/ATAC/rasqual/output/naive_100kb/naive_50kb --cis_dist 5e4 --eigenMT_path ~/software/eigenMT/eigenMT.py --external"
cat ../macrophage-gxe-study/macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir ../macrophage-gxe-study/results/SL1344/eigenMT/input/ --genepos results/ATAC/eigenMT/input/gene_positions.txt --QTL results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.eigenMT_input.txt --out_prefix results/ATAC/rasqual/output/IFNg_100kb/IFNg_50kb --cis_dist 5e4 --eigenMT_path ~/software/eigenMT/eigenMT.py --external"
cat ../macrophage-gxe-study/macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir ../macrophage-gxe-study/results/SL1344/eigenMT/input/ --genepos results/ATAC/eigenMT/input/gene_positions.txt --QTL results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.eigenMT_input.txt --out_prefix results/ATAC/rasqual/output/SL1344_100kb/SL1344_50kb --cis_dist 5e4 --eigenMT_path ~/software/eigenMT/eigenMT.py --external"
cat ../macrophage-gxe-study/macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir ../macrophage-gxe-study/results/SL1344/eigenMT/input/ --genepos results/ATAC/eigenMT/input/gene_positions.txt --QTL results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.eigenMT_input.txt --out_prefix results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_50kb --cis_dist 5e4 --eigenMT_path ~/software/eigenMT/eigenMT.py --external"

#Merge eigenMT output 
cat results/ATAC/rasqual/output/naive_100kb/naive_50kb.chr_*.eigenMT.txt | grep -v snps > results/ATAC/rasqual/output/naive_100kb/naive_50kb.eigenMT.txt
cat results/ATAC/rasqual/output/IFNg_100kb/IFNg_50kb.chr_*.eigenMT.txt | grep -v snps > results/ATAC/rasqual/output/IFNg_100kb/IFNg_50kb.eigenMT.txt
cat results/ATAC/rasqual/output/SL1344_100kb/SL1344_50kb.chr_*.eigenMT.txt | grep -v snps > results/ATAC/rasqual/output/SL1344_100kb/SL1344_50kb.eigenMT.txt
cat results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_50kb.chr_*.eigenMT.txt | grep -v snps > results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_50kb.eigenMT.txt

#Index text files with tabix
tabix -s3 -b4 -e4 -f results/ATAC/rasqual/output_fixed_gt/naive_100kb/naive_100kb.sorted.txt.gz &
tabix -s3 -b4 -e4 -f results/ATAC/rasqual/output_fixed_gt/IFNg_100kb/IFNg_100kb.sorted.txt.gz &
tabix -s3 -b4 -e4 -f results/ATAC/rasqual/output_fixed_gt/SL1344_100kb/SL1344_100kb.sorted.txt.gz &
tabix -s3 -b4 -e4 -f results/ATAC/rasqual/output_fixed_gt/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz &




