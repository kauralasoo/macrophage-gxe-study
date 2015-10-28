#Extract relevant genotypes from the vcf file
bcftools view -O z -S results/SL1344/fastqtl/input/genotype_list.txt genotypes/SL1344/imputed_20151005/imputed.59_samples.snps_indels.INFO_08.vcf.gz > results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz
#Add unique names to unnamed SNPs
bcftools annotate -O z --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz > results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz
#Remove two SNPs with duplicate records
bcftools view -O z -e 'ID=@genotypes/SL1344/imputed_20151005/duplicate_snps.txt' results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz > results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.vcf.gz
#Index the vcf file
tabix -p vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.vcf.gz 

#Compress and index expression files
bgzip results/SL1344/fastqtl/input/naive.expression.txt && tabix -p bed results/SL1344/fastqtl/input/naive.expression.txt.gz
bgzip results/SL1344/fastqtl/input/IFNg.expression.txt && tabix -p bed results/SL1344/fastqtl/input/IFNg.expression.txt.gz
bgzip results/SL1344/fastqtl/input/SL1344.expression.txt && tabix -p bed results/SL1344/fastqtl/input/SL1344.expression.txt.gz
bgzip results/SL1344/fastqtl/input/IFNg_SL1344.expression.txt && tabix -p bed results/SL1344/fastqtl/input/IFNg_SL1344.expression.txt.gz

#Run FastQTL on each condition
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.vcf.gz --bed results/SL1344/fastqtl/input/naive.expression.txt.gz --cov results/SL1344/fastqtl/input/naive.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/naive_perm --execute True"
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.vcf.gz --bed results/SL1344/fastqtl/input/IFNg.expression.txt.gz --cov results/SL1344/fastqtl/input/IFNg.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/IFNg_perm --execute True"
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.vcf.gz --bed results/SL1344/fastqtl/input/SL1344.expression.txt.gz --cov results/SL1344/fastqtl/input/SL1344.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/SL1344_perm --execute True"
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.no_dup.vcf.gz --bed results/SL1344/fastqtl/input/IFNg_SL1344.expression.txt.gz --cov results/SL1344/fastqtl/input/IFNg_SL1344.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/IFNg_SL1344_perm --execute True"

#Merge chunks into single files
zcat results/SL1344/fastqtl/output/naive_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/naive_permuted.txt.gz
zcat results/SL1344/fastqtl/output/IFNg_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_permuted.txt.gz
zcat results/SL1344/fastqtl/output/SL1344_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/SL1344_permuted.txt.gz
zcat results/SL1344/fastqtl/output/IFNg_SL1344_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_SL1344_permuted.txt.gz

