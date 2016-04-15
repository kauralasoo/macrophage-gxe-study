#Compress and index expression files
bgzip results/SL1344/fastqtl/input/naive.expression_cqn.txt && tabix -p bed results/SL1344/fastqtl/input/naive.expression_cqn.txt.gz
bgzip results/SL1344/fastqtl/input/IFNg.expression_cqn.txt && tabix -p bed results/SL1344/fastqtl/input/IFNg.expression_cqn.txt.gz
bgzip results/SL1344/fastqtl/input/SL1344.expression_cqn.txt && tabix -p bed results/SL1344/fastqtl/input/SL1344.expression_cqn.txt.gz
bgzip results/SL1344/fastqtl/input/IFNg_SL1344.expression_cqn.txt && tabix -p bed results/SL1344/fastqtl/input/IFNg_SL1344.expression_cqn.txt.gz

#Run FastQTL on each condition
cat results/SL1344/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --bed results/SL1344/fastqtl/input/naive.expression_cqn.txt.gz --cov results/SL1344/fastqtl/input/naive.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/naive_500kb_cqn_perm --execute True"
cat results/SL1344/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --bed results/SL1344/fastqtl/input/IFNg.expression_cqn.txt.gz --cov results/SL1344/fastqtl/input/IFNg.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/IFNg_500kb_cqn_perm --execute True"
cat results/SL1344/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --bed results/SL1344/fastqtl/input/SL1344.expression_cqn.txt.gz --cov results/SL1344/fastqtl/input/SL1344.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/SL1344_500kb_cqn_perm --execute True"
cat results/SL1344/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --bed results/SL1344/fastqtl/input/IFNg_SL1344.expression_cqn.txt.gz --cov results/SL1344/fastqtl/input/IFNg_SL1344.covariates.txt --W 500000 --permute '100 10000' --out results/SL1344/fastqtl/output/IFNg_SL1344_500kb_cqn_perm --execute True"


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
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --fastqtl results/SL1344/fastqtl/output/naive_pvalues.txt.gz | bgzip > results/SL1344/fastqtl/output/naive_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --fastqtl results/SL1344/fastqtl/output/IFNg_pvalues.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --fastqtl results/SL1344/fastqtl/output/SL1344_pvalues.txt.gz | bgzip > results/SL1344/fastqtl/output/SL1344_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/fastqtl/input/fastqtl_genotypes.INFO_08.named.vcf.gz --fastqtl results/SL1344/fastqtl/output/IFNg_SL1344_pvalues.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_SL1344_pvalues.coords.txt.gz"

#Extract p-value for a subset of genes ENSG00000185245 ENSG00000079785 ENSG00000144228 ENSG00000173992
echo "CHR BP SNP P\n" && zgrep ENSG00000109861 results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " > results/SL1344/eQTLs/example_loci/CTSC/CTSC_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/CCDC9_pvalues.gwas && zgrep ENSG00000105321 results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/CCDC9_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/ITGB3BP_pvalues.gwas && zgrep ENSG00000142856 results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/ITGB3BP_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/FCGR3A_pvalues.gwas && zgrep ENSG00000203747 results/SL1344/fastqtl/output/naive_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/FCGR3A_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/FCGR3B_pvalues.gwas && zgrep ENSG00000162747 results/SL1344/fastqtl/output/naive_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/FCGR3B_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/GOLGA7_pvalues.gwas && zgrep ENSG00000147533 results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/GOLGA7_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/GP1BA_pvalues.gwas && zgrep ENSG00000185245 results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/GP1BA_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/DDX11_pvalues.gwas && zgrep ENSG00000013573 results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/DDX11_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/SPOPL_pvalues.gwas && zgrep ENSG00000144228 results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/SPOPL_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/CCS_pvalues.gwas && zgrep ENSG00000173992 results/SL1344/fastqtl/output/naive_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/CCS_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/RGS14_pvalues.gwas && zgrep ENSG00000169220 results/SL1344/fastqtl/output/IFNg_SL1344_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/RGS14_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/TMEM229B_pvalues.gwas && zgrep ENSG00000198133 results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/TMEM229B_pvalues.gwas
printf "CHR BP SNP P\n" > results/SL1344/eQTLs/example_loci/TLR1_pvalues.gwas && zgrep ENSG00000174125 results/SL1344/fastqtl/output/SL1344_pvalues.coords.txt.gz | cut -f2,3,4,6 -d " " >> results/SL1344/eQTLs/example_loci/TLR1_pvalues.gwas


python ~/software/utils/fastqtl/fastqtlExtractGenePvals.py --fastqtl results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz --gene_id ENSG00000144228 > SPOPL.gwas

#IL2RA pvalues
python ~/software/utils/fastqtl/fastqtlExtractGenePvals.py --fastqtl results/SL1344/fastqtl/output/IFNg_pvalues.coords.txt.gz --gene_id ENSG00000134460 > results/SL1344/eQTLs/example_loci/IL2RA_IFNg_pvalues.gwas
python ~/software/utils/fastqtl/fastqtlExtractGenePvals.py --fastqtl results/SL1344/fastqtl/output/naive_pvalues.coords.txt.gz --gene_id ENSG00000134460 > results/SL1344/eQTLs/example_loci/IL2RA_naive_pvalues.gwas
python ~/software/utils/fastqtl/fastqtlExtractGenePvals.py --fastqtl results/SL1344/fastqtl/output/SL1344_pvalues.coords.txt.gz --gene_id ENSG00000134460 > results/SL1344/eQTLs/example_loci/IL2RA_SL1344_pvalues.gwas
python ~/software/utils/fastqtl/fastqtlExtractGenePvals.py --fastqtl results/SL1344/fastqtl/output/IFNg_SL1344_pvalues.coords.txt.gz --gene_id ENSG00000134460 > results/SL1344/eQTLs/example_loci/IL2RA_IFNg_SL1344_pvalues.gwas



