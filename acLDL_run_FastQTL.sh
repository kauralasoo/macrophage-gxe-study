#Extract relevant genotypes from the vcf file
bcftools view -S results/acLDL/fastqtl/input/genotype_list.txt genotypes/acLDL/imputed_20151005/imputed.45_samples.snps_indels.INFO_08.vcf.gz | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' - | bcftools view -O z -e 'ID=@genotypes/SL1344/imputed_20151005/duplicate_snps.txt' -  > results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz
#Index the vcf file
tabix -p vcf results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz 

#Compress and index expression files
bgzip results/acLDL/fastqtl/input/Ctrl.expression.txt && tabix -p bed results/acLDL/fastqtl/input/Ctrl.expression.txt.gz
bgzip results/acLDL/fastqtl/input/AcLDL.expression.txt && tabix -p bed results/acLDL/fastqtl/input/AcLDL.expression.txt.gz

#Run FastQTL on each condition
cat results/acLDL/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/runFastQTL.py --vcf results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz --bed results/acLDL/fastqtl/input/Ctrl.expression.txt.gz --cov results/acLDL/fastqtl/input/Ctrl.covariates.txt --W 500000 --permute '100 10000' --out results/acLDL/fastqtl/output/Ctrl_perm --execute True"
cat results/acLDL/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/runFastQTL.py --vcf results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz --bed results/acLDL/fastqtl/input/AcLDL.expression.txt.gz --cov results/acLDL/fastqtl/input/AcLDL.covariates.txt --W 500000 --permute '100 10000' --out results/acLDL/fastqtl/output/AcLDL_perm --execute True"

#Merge chunks into single files
zcat results/acLDL/fastqtl/output/Ctrl_perm.chunk_*.txt.gz | bgzip > results/acLDL/fastqtl/output/Ctrl_permuted.txt.gz
zcat results/acLDL/fastqtl/output/AcLDL_perm.chunk_*.txt.gz | bgzip > results/acLDL/fastqtl/output/AcLDL_permuted.txt.gz

#Run on all SNPs 
cat results/acLDL/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/runFastQTL.py --vcf results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz --bed results/acLDL/fastqtl/input/Ctrl.expression.txt.gz --cov results/acLDL/fastqtl/input/Ctrl.covariates.txt --W 500000 --out results/acLDL/fastqtl/output/Ctrl_full --execute True"
cat results/acLDL/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/runFastQTL.py --vcf results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz --bed results/acLDL/fastqtl/input/AcLDL.expression.txt.gz --cov results/acLDL/fastqtl/input/AcLDL.covariates.txt --W 500000 --out results/acLDL/fastqtl/output/AcLDL_full --execute True"

#Merge chunks into a single file
zcat results/acLDL/fastqtl/output/Ctrl_full.chunk_*.txt.gz | bgzip > results/acLDL/fastqtl/output/Ctrl_pvalues.txt.gz
zcat results/acLDL/fastqtl/output/AcLDL_full.chunk_*.txt.gz | bgzip > results/acLDL/fastqtl/output/AcLDL_pvalues.txt.gz

#Remove the chunks
rm results/acLDL/fastqtl/output/Ctrl_full.chunk_*.txt.gz 
rm results/acLDL/fastqtl/output/AcLDL_full.chunk_*.txt.gz 