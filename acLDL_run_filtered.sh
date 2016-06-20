#Extract genotype ids for each condition
cut -f2 results/acLDL/rasqual/input_filtered/Ctrl.sg_map.txt > results/acLDL/rasqual/input_filtered/Ctrl.genotypes.txt
cut -f2 results/acLDL/rasqual/input_filtered/AcLDL.sg_map.txt > results/acLDL/rasqual/input_filtered/AcLDL.genotypes.txt

#Extract samples from the global VCF file
bcftools view -Oz -S results/acLDL/rasqual/input_filtered/AcLDL.genotypes.txt genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/acLDL/rasqual/input_filtered/AcLDL.vcf &
bcftools view -Oz -S results/acLDL/rasqual/input_filtered/Ctrl.genotypes.txt genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/acLDL/rasqual/input_filtered/Ctrl.vcf &

#Add ASE counts into the VCF file
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 32000 --jobname vcfAddASE --queue hugemem --command "python ~/software/rasqual/scripts/vcfAddASE.py --ASEcounts results/acLDL/combined_ASE_counts.txt --ASESampleGenotypeMap results/acLDL/rasqual/input_filtered/Ctrl.sg_map.txt --VCFfile results/acLDL/rasqual/input_filtered/Ctrl.vcf | bgzip > results/acLDL/rasqual/input_filtered/Ctrl.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 32000 --jobname vcfAddASE --queue hugemem --command "python ~/software/rasqual/scripts/vcfAddASE.py --ASEcounts results/acLDL/combined_ASE_counts.txt --ASESampleGenotypeMap results/acLDL/rasqual/input_filtered/AcLDL.sg_map.txt --VCFfile results/acLDL/rasqual/input_filtered/AcLDL.vcf | bgzip > results/acLDL/rasqual/input_filtered/AcLDL.ASE.vcf.gz"

#Index VCF files
tabix -p vcf results/acLDL/rasqual/input_filtered/Ctrl.ASE.vcf.gz 
tabix -p vcf results/acLDL/rasqual/input_filtered/AcLDL.ASE.vcf.gz 

#RASQUAL (2 PEER covariates + sex, library sizes + GC), 500kb
cat results/acLDL/rasqual/input_filtered/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_Ctrl --ncores 8 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/acLDL/rasqual/input_filtered/Ctrl.expression.bin  --offsets results/acLDL/rasqual/input_filtered/Ctrl.gc_library_size.bin --n 58 --geneids results/acLDL/rasqual/input_filtered/feature_names.txt --vcf results/acLDL/rasqual/input_filtered/Ctrl.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input_filtered/gene_snp_count_500kb.txt --outprefix results/acLDL/rasqual/output_filtered/Ctrl_500kb/batches/Ctrl_500kb --covariates results/acLDL/rasqual/input_filtered/Ctrl.covariates.bin --rasqualBin rasqual --parameters '\--force --n-threads 8' --execute True"

cat results/acLDL/rasqual/input_filtered/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_AcLDL --ncores 8 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/acLDL/rasqual/input_filtered/AcLDL.expression.bin  --offsets results/acLDL/rasqual/input_filtered/AcLDL.gc_library_size.bin --n 58 --geneids results/acLDL/rasqual/input_filtered/feature_names.txt --vcf results/acLDL/rasqual/input_filtered/AcLDL.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input_filtered/gene_snp_count_500kb.txt --outprefix results/acLDL/rasqual/output_filtered/AcLDL_500kb/batches/AcLDL_500kb --covariates results/acLDL/rasqual/input_filtered/AcLDL.covariates.bin --rasqualBin rasqual --parameters '\--force --n-threads 8' --execute True"

#Merge batches
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/acLDL/rasqual/output_filtered/Ctrl_500kb/batches/Ctrl_500kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/acLDL/rasqual/output_filtered/AcLDL_500kb/batches/AcLDL_500kb"

#Move merged file to parent dir
mv results/acLDL/rasqual/output_filtered/Ctrl_500kb/batches/Ctrl_500kb.txt results/acLDL/rasqual/output_filtered/Ctrl_500kb/
mv results/acLDL/rasqual/output_filtered/AcLDL_500kb/batches/AcLDL_500kb.txt results/acLDL/rasqual/output_filtered/AcLDL_500kb/

#Identify completed genes
cut -f1 results/acLDL/rasqual/output_filtered/Ctrl_500kb/Ctrl_500kb.txt | uniq > results/acLDL/rasqual/output_filtered/Ctrl_500kb/Ctrl_500kb.completed_ids.txt
cut -f1 results/acLDL/rasqual/output_filtered/AcLDL_500kb/AcLDL_500kb.txt | uniq > results/acLDL/rasqual/output_filtered/AcLDL_500kb/AcLDL_500kb.completed_ids.txt

# Sort and filter merged p-values
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED results/acLDL/rasqual/output_filtered/Ctrl_500kb/Ctrl_500kb.txt | sort -k3,3 -k4,4n | bgzip > results/acLDL/rasqual/output_filtered/Ctrl_500kb/Ctrl_500kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED results/acLDL/rasqual/output_filtered/AcLDL_500kb/AcLDL_500kb.txt | sort -k3,3 -k4,4n | bgzip > results/acLDL/rasqual/output_filtered/AcLDL_500kb/AcLDL_500kb.sorted.txt.gz"

#Index the output files using Tabix
tabix -s3 -b4 -e4 -f results/acLDL/rasqual/output_filtered/AcLDL_500kb/AcLDL_500kb.sorted.txt.gz
tabix -s3 -b4 -e4 -f results/acLDL/rasqual/output_filtered/Ctrl_500kb/Ctrl_500kb.sorted.txt.gz


### eigenMT ####
#Convert rasqual output into format suitable for eigenMT
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/acLDL/rasqual/output_filtered/Ctrl_500kb/Ctrl_500kb.txt > results/acLDL/rasqual/output_filtered/Ctrl_500kb/Ctrl_500kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/acLDL/rasqual/output_filtered/AcLDL_500kb/AcLDL_500kb.txt > results/acLDL/rasqual/output_filtered/AcLDL_500kb/AcLDL_500kb.eigenMT_input.txt"

#Run eigenMT on on each chromosome
cat macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 4 --jobname eigenMTbyChromosome2 --command "python ~/software/utils/eigenMTbyChromosome.py --genepos results/acLDL/rasqual/input/gene_positions.txt --chromosome_dir genotypes/acLDL/imputed_20151005/chromosomes_INFO_07/ --QTL results/acLDL/rasqual/output_filtered/Ctrl_500kb/Ctrl_500kb.eigenMT_input.txt --out_prefix results/acLDL/rasqual/output_filtered/Ctrl_500kb/Ctrl_500kb --cis_dist 1e7 --eigenMT_path ~/software/eigenMT/eigenMT.py"
cat macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 4 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --genepos results/acLDL/rasqual/input/gene_positions.txt --chromosome_dir genotypes/acLDL/imputed_20151005/chromosomes_INFO_07/ --QTL results/acLDL/rasqual/output_filtered/AcLDL_500kb/AcLDL_500kb.eigenMT_input.txt --out_prefix results/acLDL/rasqual/output_filtered/AcLDL_500kb/AcLDL_500kb --cis_dist 1e7 --eigenMT_path ~/software/eigenMT/eigenMT.py"



