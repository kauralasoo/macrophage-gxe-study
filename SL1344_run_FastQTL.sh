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
zcat results/SL1344/fastqtl/output/naive_500kb_cqn_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/naive_500kb_permuted.txt.gz
zcat results/SL1344/fastqtl/output/IFNg_500kb_cqn_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_500kb_permuted.txt.gz
zcat results/SL1344/fastqtl/output/SL1344_500kb_cqn_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/SL1344_500kb_permuted.txt.gz
zcat results/SL1344/fastqtl/output/IFNg_SL1344_500kb_cqn_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_SL1344_500kb_permuted.txt.gz

#Remove chunks
rm results/SL1344/fastqtl/output/*.chunk_*


#Get full p-values from fastQTL
cat results/SL1344/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --bed results/SL1344/fastqtl/input/naive.expression_cqn.txt.gz --cov results/SL1344/fastqtl/input/naive.covariates.txt --W 500000 --out results/SL1344/fastqtl/output/naive_500kb_full --execute True"
cat results/SL1344/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --bed results/SL1344/fastqtl/input/IFNg.expression_cqn.txt.gz --cov results/SL1344/fastqtl/input/IFNg.covariates.txt --W 500000 --out results/SL1344/fastqtl/output/IFNg_500kb_full --execute True"
cat results/SL1344/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --bed results/SL1344/fastqtl/input/SL1344.expression_cqn.txt.gz --cov results/SL1344/fastqtl/input/SL1344.covariates.txt --W 500000 --out results/SL1344/fastqtl/output/SL1344_500kb_full --execute True"
cat results/SL1344/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --bed results/SL1344/fastqtl/input/IFNg_SL1344.expression_cqn.txt.gz --cov results/SL1344/fastqtl/input/IFNg_SL1344.covariates.txt --W 500000 --out results/SL1344/fastqtl/output/IFNg_SL1344_500kb_full --execute True"


#Merge chunks into single files
zcat results/SL1344/fastqtl/output/naive_500kb_full.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/naive_500kb_pvalues.txt.gz
zcat results/SL1344/fastqtl/output/IFNg_500kb_full.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_500kb_pvalues.txt.gz
zcat results/SL1344/fastqtl/output/SL1344_500kb_full.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/SL1344_500kb_pvalues.txt.gz
zcat results/SL1344/fastqtl/output/IFNg_SL1344_500kb_full.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_SL1344_500kb_pvalues.txt.gz

#Remove chunks
rm results/SL1344/fastqtl/output/*.chunk_*

#Add SNP coordinates
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --fastqtl results/SL1344/fastqtl/output/naive_500kb_pvalues.txt.gz | bgzip > results/SL1344/fastqtl/output/naive_500kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --fastqtl results/SL1344/fastqtl/output/IFNg_500kb_pvalues.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_500kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --fastqtl results/SL1344/fastqtl/output/SL1344_500kb_pvalues.txt.gz | bgzip > results/SL1344/fastqtl/output/SL1344_500kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --fastqtl results/SL1344/fastqtl/output/IFNg_SL1344_500kb_pvalues.txt.gz | bgzip > results/SL1344/fastqtl/output/IFNg_SL1344_500kb_pvalues.coords.txt.gz"

#Sort files by SNP coordinates
#awk command is necessary to change field separator from space to tab
zcat results/SL1344/fastqtl/output/naive_500kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/SL1344/fastqtl/output/naive_500kb_pvalues.sorted.txt.gz &
zcat results/SL1344/fastqtl/output/IFNg_500kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/SL1344/fastqtl/output/IFNg_500kb_pvalues.sorted.txt.gz &
zcat results/SL1344/fastqtl/output/SL1344_500kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/SL1344/fastqtl/output/SL1344_500kb_pvalues.sorted.txt.gz &
zcat results/SL1344/fastqtl/output/IFNg_SL1344_500kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/SL1344/fastqtl/output/IFNg_SL1344_500kb_pvalues.sorted.txt.gz &

#Index the output files using Tabix
tabix -s2 -b3 -e3 -f results/SL1344/fastqtl/output/naive_500kb_pvalues.sorted.txt.gz
tabix -s2 -b3 -e3 -f results/SL1344/fastqtl/output/IFNg_500kb_pvalues.sorted.txt.gz
tabix -s2 -b3 -e3 -f results/SL1344/fastqtl/output/SL1344_500kb_pvalues.sorted.txt.gz
tabix -s2 -b3 -e3 -f results/SL1344/fastqtl/output/IFNg_SL1344_500kb_pvalues.sorted.txt.gz

#Remove intermediate files
rm results/SL1344/fastqtl/output/*pvalues.coords.txt.gz
rm results/SL1344/fastqtl/output/*pvalues.txt.gz


###### Map splicing QTLs using LeafCutter ######

#Construct junction files
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | head -n 100 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamToJunctions --command "python ~/software/utils/bam/bamToJunctions.py --indir STAR/SL1344/ --outdir STAR/SL1344/ --leafCutterDir ~/software/leafcutter/ --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | head -n 200 | tail -n 100 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamToJunctions --command "python ~/software/utils/bam/bamToJunctions.py --indir STAR/SL1344/ --outdir STAR/SL1344/ --leafCutterDir ~/software/leafcutter/ --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | tail -n 160 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamToJunctions --command "python ~/software/utils/bam/bamToJunctions.py --indir STAR/SL1344/ --outdir STAR/SL1344/ --leafCutterDir ~/software/leafcutter/ --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Cluster junctions
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o FarmOut/leafcutter_cluster.%J.jobout "python ~/software/leafcutter/clustering/leafcutter_cluster.py -j results/SL1344/leafcutter/junction_files.txt -r results/SL1344/leafcutter/ -m 50 -l 500000"

#Split leafcutter output into cluster counts and intron counts
python ~/software/utils/leafcutter/leafcutter_split_counts.py --leafcutter_out results/SL1344/leafcutter/leafcutter_perind.counts.gz --outprefix leafcutter



