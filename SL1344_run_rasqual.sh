#Add read group to the BAM files to make them work with GATK
cat macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamAddRG --command "python ~/software/utils/bam/bamAddRG.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --outsuffix .Aligned.sortedByCoord.RG.bam --execute True"
cat macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamAddRG --command "python ~/software/utils/bam/bamAddRG.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --outsuffix .Aligned.sortedByCoord.RG.bam --execute True"


#Index bams
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Index the vcf file
tabix -p vcf genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz
bcftools index genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz

#Use ASEReadCounter to count allele-specific expression
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | python ~/software/utils/submitJobs.py --MEM 6000 --jobname bamCountASE --command "python ~/software/rasqual/scripts/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True --Xmx 4g"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map_2.txt | python ~/software/utils/submitJobs.py --MEM 10000 --jobname bamCountASE --command "python ~/software/rasqual/scripts/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True --Xmx 8g"

#Construct a sample-sample map for meregeASECounts.py script
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | awk -v OFS='\t' '{print $1, $1}' > results/SL1344/rasqual/input/sample_sample_map.txt

#Merge all allele-specific counts into one matrix
echo "mergeASECounts" | python ~/software/utils/submitJobs.py --MEM 27000 --jobname mergeASECounts --command "python ~/software/rasqual/scripts/mergeASECounts.py --sample_list results/SL1344/rasqual/input/sample_sample_map.txt --indir STAR/SL1344 --suffix .ASEcounts | bgzip > results/SL1344/combined_ASE_counts.txt"

#Sort and index the ASE counts file
(zcat results/SL1344/combined_ASE_counts.txt.gz | head -n1 && zcat results/SL1344/combined_ASE_counts.txt.gz | tail -n+2 | sort -k1,1 -k2,2n) | bgzip > combined_ASE_counts.sorted.txt.gz
tabix -s 1 -b 2 -e 2 -S 1 combined_ASE_counts.sorted.txt.gz 

#Extract genotype ids for each condition
cut -f2 results/SL1344/rasqual/input/naive.sg_map.txt > results/SL1344/rasqual/input/naive.genotypes.txt
cut -f2 results/SL1344/rasqual/input/IFNg.sg_map.txt > results/SL1344/rasqual/input/IFNg.genotypes.txt
cut -f2 results/SL1344/rasqual/input/SL1344.sg_map.txt > results/SL1344/rasqual/input/SL1344.genotypes.txt
cut -f2 results/SL1344/rasqual/input/IFNg_SL1344.sg_map.txt > results/SL1344/rasqual/input/IFNg_SL1344.genotypes.txt

#Extract samples from the global VCF file
bcftools view -Oz -S results/SL1344/rasqual/input/naive.genotypes.txt genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/SL1344/rasqual/input/naive.vcf &
bcftools view -Oz -S results/SL1344/rasqual/input/IFNg.genotypes.txt genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/SL1344/rasqual/input/IFNg.vcf &
bcftools view -Oz -S results/SL1344/rasqual/input/SL1344.genotypes.txt genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/SL1344/rasqual/input/SL1344.vcf &
bcftools view -Oz -S results/SL1344/rasqual/input/IFNg_SL1344.genotypes.txt genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/SL1344/rasqual/input/IFNg_SL1344.vcf &

#Add ASE counts into the VCF file
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 42000 --jobname vcfAddASE --queue hugemem --command "python ~/software/rasqual/scripts/vcfAddASE.py --ASEcounts results/SL1344/combined_ASE_counts.txt --ASESampleGenotypeMap results/SL1344/rasqual/input/naive.sg_map.txt --VCFfile results/SL1344/rasqual/input/naive.vcf | bgzip > results/SL1344/rasqual/input/naive.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 42000 --jobname vcfAddASE --queue hugemem --command "python ~/software/rasqual/scripts/vcfAddASE.py --ASEcounts results/SL1344/combined_ASE_counts.txt --ASESampleGenotypeMap results/SL1344/rasqual/input/IFNg.sg_map.txt --VCFfile results/SL1344/rasqual/input/IFNg.vcf | bgzip > results/SL1344/rasqual/input/IFNg.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 42000 --jobname vcfAddASE --queue hugemem --command "python ~/software/rasqual/scripts/vcfAddASE.py --ASEcounts results/SL1344/combined_ASE_counts.txt --ASESampleGenotypeMap results/SL1344/rasqual/input/SL1344.sg_map.txt --VCFfile results/SL1344/rasqual/input/SL1344.vcf | bgzip > results/SL1344/rasqual/input/SL1344.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 42000 --jobname vcfAddASE --queue hugemem --command "python ~/software/rasqual/scripts/vcfAddASE.py --ASEcounts results/SL1344/combined_ASE_counts.txt --ASESampleGenotypeMap results/SL1344/rasqual/input/IFNg_SL1344.sg_map.txt --VCFfile results/SL1344/rasqual/input/IFNg_SL1344.vcf | bgzip > results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz"

#Index VCF files using tabix
tabix -p vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz &
tabix -p vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz &
tabix -p vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz &
tabix -p vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz &

#RUN rasqual on all peaks
#naive
cat results/SL1344/rasqual/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 8 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 84 --geneids results/SL1344/rasqual/input/gene_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/naive_500kb/batches/naive_500kb --covariates results/SL1344/rasqual/input/naive.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 8'"
cat results/SL1344/rasqual/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_IFNg --ncores 8 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/SL1344/rasqual/input/IFNg.expression.bin --offsets results/SL1344/rasqual/input/IFNg.gc_library_size.bin --n 84 --geneids results/SL1344/rasqual/input/gene_names.txt --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/IFNg_500kb/batches/IFNg_500kb --covariates results/SL1344/rasqual/input/IFNg.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 8'"
cat results/SL1344/rasqual/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_SL1344 --ncores 8 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/SL1344/rasqual/input/SL1344.expression.bin --offsets results/SL1344/rasqual/input/SL1344.gc_library_size.bin --n 84 --geneids results/SL1344/rasqual/input/gene_names.txt --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/SL1344_500kb/batches/SL1344_500kb --covariates results/SL1344/rasqual/input/SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 8'"
cat results/SL1344/rasqual/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_IFNg_SL1344 --ncores 8 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/SL1344/rasqual/input/IFNg_SL1344.expression.bin --offsets results/SL1344/rasqual/input/IFNg_SL1344.gc_library_size.bin --n 84 --geneids results/SL1344/rasqual/input/gene_names.txt --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/IFNg_SL1344_500kb/batches/IFNg_SL1344_500kb --covariates results/SL1344/rasqual/input/IFNg_SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 8'"

#Rerun failed batches
cat results/SL1344/rasqual/input/naive_failed_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --queue basement --ncores 8 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 84 --geneids results/SL1344/rasqual/input/gene_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/naive_500kb/batches/naive_500kb --covariates results/SL1344/rasqual/input/naive.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 8 --imputation-quality-fsnp 0.9'"
cat results/SL1344/rasqual/input/IFNg_failed_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_IFNg --queue basement --ncores 8 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/SL1344/rasqual/input/IFNg.expression.bin --offsets results/SL1344/rasqual/input/IFNg.gc_library_size.bin --n 84 --geneids results/SL1344/rasqual/input/gene_names.txt --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/IFNg_500kb/batches/IFNg_500kb --covariates results/SL1344/rasqual/input/IFNg.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 8 --imputation-quality-fsnp 0.9'"
cat results/SL1344/rasqual/input/SL1344_failed_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_SL1344 --queue basement --ncores 8 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/SL1344/rasqual/input/SL1344.expression.bin --offsets results/SL1344/rasqual/input/SL1344.gc_library_size.bin --n 84 --geneids results/SL1344/rasqual/input/gene_names.txt --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/SL1344_500kb/batches/SL1344_500kb --covariates results/SL1344/rasqual/input/SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 8 --imputation-quality-fsnp 0.9'"
cat results/SL1344/rasqual/input/IFNg_SL1344_failed_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_IFNg_SL1344 --queue basement --ncores 8 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/SL1344/rasqual/input/IFNg_SL1344.expression.bin --offsets results/SL1344/rasqual/input/IFNg_SL1344.gc_library_size.bin --n 84 --geneids results/SL1344/rasqual/input/gene_names.txt --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/IFNg_SL1344_500kb/batches/IFNg_SL1344_500kb --covariates results/SL1344/rasqual/input/IFNg_SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 8 --imputation-quality-fsnp 0.9'"


#Merge batches into single files
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/SL1344/rasqual/output/naive_500kb/batches/naive_500kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/SL1344/rasqual/output/IFNg_500kb/batches/IFNg_500kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/SL1344/rasqual/output/SL1344_500kb/batches/SL1344_500kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/SL1344/rasqual/output/IFNg_SL1344_500kb/batches/IFNg_SL1344_500kb"

#Move merged files to parent directory
mv results/SL1344/rasqual/output/naive_500kb/batches/naive_500kb.txt results/SL1344/rasqual/output/naive_500kb/
mv results/SL1344/rasqual/output/IFNg_500kb/batches/IFNg_500kb.txt results/SL1344/rasqual/output/IFNg_500kb/
mv results/SL1344/rasqual/output/SL1344_500kb/batches/SL1344_500kb.txt results/SL1344/rasqual/output/SL1344_500kb/
mv results/SL1344/rasqual/output/IFNg_SL1344_500kb/batches/IFNg_SL1344_500kb.txt results/SL1344/rasqual/output/IFNg_SL1344_500kb/

#Extract completed gene ids
cut -f1 results/SL1344/rasqual/output/naive_500kb/naive_500kb.txt | uniq > results/SL1344/rasqual/output/naive_500kb/naive_500kb.completed_ids.txt &
cut -f1 results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.txt | uniq > results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.completed_ids.txt &
cut -f1 results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.txt | uniq > results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.completed_ids.txt &
cut -f1 results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.txt | uniq > results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.completed_ids.txt &

# Sort and filter merged p-values
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED results/SL1344/rasqual/output/naive_500kb/naive_500kb.txt | sort -k3,3 -k4,4n | bgzip > results/SL1344/rasqual/output/naive_500kb/naive_500kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.txt | sort -k3,3 -k4,4n | bgzip > results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.txt | sort -k3,3 -k4,4n | bgzip > results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.txt | sort -k3,3 -k4,4n | bgzip > results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.sorted.txt.gz"

#Index the output files using Tabix
tabix -s3 -b4 -e4 -f results/SL1344/rasqual/output/naive_500kb/naive_500kb.sorted.txt.gz
tabix -s3 -b4 -e4 -f results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.sorted.txt.gz
tabix -s3 -b4 -e4 -f results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.sorted.txt.gz
tabix -s3 -b4 -e4 -f results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.sorted.txt.gz

### eigenMT ####
#Split vcf into chromosomes
cat ../../../macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/vcf/vcfSplitByChromosome.py --vcf imputed.86_samples.sorted.filtered.named.INFO_07.vcf.gz --outdir chromosomes_INFO_07/ --execute False

#Convert vcfs to GDS
/software/R-3.1.2/bin/Rscript ~/software/utils/vcf/vcfToGds.R --vcf-directory chromosomes_IFNO_07 --chr-list ../../../macrophage-gxe-study/data/sample_lists/chromosome_list.txt

#Convert rasqual output into format suitable for eigenMT
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/SL1344/rasqual/output/naive_500kb/naive_500kb.txt > results/SL1344/rasqual/output/naive_500kb/naive_500kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.txt > results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.txt > results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.txt > results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.eigenMT_input.txt"

#Run eigenMT chromosome-by-chromosme
cat macrophage-gxe-study/data/sample_lists/chromosome_list.txt |  python ~/software/utils/submitJobs.py --MEM 2000 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir results/SL1344/eigenMT/input/ --genepos results/SL1344/eigenMT/input/gene_positions.txt --QTL results/SL1344/rasqual/output/naive_500kb/naive_500kb.eigenMT_input.txt --out_prefix results/SL1344/rasqual/output/naive_500kb/naive_500kb --cis_dist 5e7 --eigenMT_path ~/software/eigenMT/eigenMT.py"
cat macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir results/SL1344/eigenMT/input/ --genepos results/SL1344/eigenMT/input/gene_positions.txt --QTL results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.eigenMT_input.txt --out_prefix results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb --cis_dist 5e7  --eigenMT_path ~/software/eigenMT/eigenMT.py"
cat macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir results/SL1344/eigenMT/input/ --genepos results/SL1344/eigenMT/input/gene_positions.txt --QTL results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.eigenMT_input.txt --out_prefix results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb --cis_dist 5e7  --eigenMT_path ~/software/eigenMT/eigenMT.py"
cat macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir results/SL1344/eigenMT/input/ --genepos results/SL1344/eigenMT/input/gene_positions.txt --QTL results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.eigenMT_input.txt --out_prefix results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb --cis_dist 5e7  --eigenMT_path ~/software/eigenMT/eigenMT.py"

#Concat all eigenMT outputs
cat results/SL1344/rasqual/output/naive_500kb/naive_500kb.chr_*.eigenMT.txt | grep -v snps > results/SL1344/rasqual/output/naive_500kb/naive_500kb.eigenMT.txt
cat results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.chr_*.eigenMT.txt | grep -v snps > results/SL1344/rasqual/output/IFNg_500kb/IFNg_500kb.eigenMT.txt
cat results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.chr_*.eigenMT.txt | grep -v snps > results/SL1344/rasqual/output/SL1344_500kb/SL1344_500kb.eigenMT.txt
cat results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.chr_*.eigenMT.txt | grep -v snps > results/SL1344/rasqual/output/IFNg_SL1344_500kb/IFNg_SL1344_500kb.eigenMT.txt



