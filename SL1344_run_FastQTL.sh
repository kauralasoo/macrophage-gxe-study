#Add unique names to unnamed SNPs and remove duplicate snps
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' genotypes/SL1344/imputed_20151005/imputed.69_samples.snps_indels.INFO_08.vcf.gz | bcftools view -O z -e 'ID=@genotypes/SL1344/imputed_20151005_old/duplicate_snps.txt' - > results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz

#Index the vcf file
tabix -p vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz 

#Compress and index expression files
bgzip results/SL1344/fastqtl/input/naive.expression.txt && tabix -p bed results/SL1344/fastqtl/input/naive.expression.txt.gz
bgzip results/SL1344/fastqtl/input/IFNg.expression.txt && tabix -p bed results/SL1344/fastqtl/input/IFNg.expression.txt.gz
bgzip results/SL1344/fastqtl/input/SL1344.expression.txt && tabix -p bed results/SL1344/fastqtl/input/SL1344.expression.txt.gz
bgzip results/SL1344/fastqtl/input/IFNg_SL1344.expression.txt && tabix -p bed results/SL1344/fastqtl/input/IFNg_SL1344.expression.txt.gz

#Run FastQTL on each condition
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --bed results/SL1344/fastqtl/input/naive.expression.txt.gz --cov results/SL1344/fastqtl/input/naive.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/naive_perm --execute True"
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --bed results/SL1344/fastqtl/input/IFNg.expression.txt.gz --cov results/SL1344/fastqtl/input/IFNg.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/IFNg_perm --execute True"
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --bed results/SL1344/fastqtl/input/SL1344.expression.txt.gz --cov results/SL1344/fastqtl/input/SL1344.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/SL1344_perm --execute True"
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --bed results/SL1344/fastqtl/input/IFNg_SL1344.expression.txt.gz --cov results/SL1344/fastqtl/input/IFNg_SL1344.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/IFNg_SL1344_perm --execute True"


#Merge chunks into single files
zcat results/SL1344/fastqtl/output/naive_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/naive_permuted.txt.gz
zcat results/SL1344/fastqtl/output/IFNg_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_permuted.txt.gz
zcat results/SL1344/fastqtl/output/SL1344_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/SL1344_permuted.txt.gz
zcat results/SL1344/fastqtl/output/IFNg_SL1344_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_SL1344_permuted.txt.gz

#Remove chunks
rm results/SL1344/fastqtl/output/*.chunk_*

#Get full p-values from fastQTL
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --bed results/SL1344/fastqtl/input/naive.expression.txt.gz --cov results/SL1344/fastqtl/input/naive.covariates.txt --W 500000 --out results/SL1344/fastqtl/output/naive_full --execute True"
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --bed results/SL1344/fastqtl/input/IFNg.expression.txt.gz --cov results/SL1344/fastqtl/input/IFNg.covariates.txt --W 500000 --out results/SL1344/fastqtl/output/IFNg_full --execute True"
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --bed results/SL1344/fastqtl/input/SL1344.expression.txt.gz --cov results/SL1344/fastqtl/input/SL1344.covariates.txt --W 500000 --out results/SL1344/fastqtl/output/SL1344_full --execute True"
cat results/SL1344/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --bed results/SL1344/fastqtl/input/IFNg_SL1344.expression.txt.gz --cov results/SL1344/fastqtl/input/IFNg_SL1344.covariates.txt --W 500000 --out results/SL1344/fastqtl/output/IFNg_SL1344_full --execute True"


#Merge chunks into single files
zcat results/SL1344/fastqtl/output/naive_full.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/naive_pvalues.txt.gz
zcat results/SL1344/fastqtl/output/IFNg_full.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_pvalues.txt.gz
zcat results/SL1344/fastqtl/output/SL1344_full.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/SL1344_pvalues.txt.gz
zcat results/SL1344/fastqtl/output/IFNg_SL1344_full.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_SL1344_pvalues.txt.gz

#Remove chunks
rm results/SL1344/fastqtl/output/*.chunk_*

#Add SNP coordinates
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --fastqtl results/SL1344/fastqtl/output/naive_pvalues.txt.gz | bgzip > results/acLDL/fastqtl/output/naive_pvalues.coords.txt.gz"


