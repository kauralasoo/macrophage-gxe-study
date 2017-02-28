#Extract relevant genotypes from the vcf file
bcftools view -S results/acLDL/fastqtl/input/genotype_list.txt genotypes/acLDL/imputed_20151005/imputed.57_samples.snps_indels.INFO_08.vcf.gz | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' - | bcftools view -O z -e 'ID=@genotypes/SL1344/imputed_20151005/duplicate_snps.txt' -  > results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz
#Index the vcf file
tabix -p vcf results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_08.vcf.gz 

#Compress and index expression files
bgzip results/acLDL/fastqtl/input/Ctrl.cqn_expression.txt && tabix -p bed results/acLDL/fastqtl/input/Ctrl.cqn_expression.txt.gz
bgzip results/acLDL/fastqtl/input/AcLDL.cqn_expression.txt && tabix -p bed results/acLDL/fastqtl/input/AcLDL.cqn_expression.txt.gz

#Run FastQTL on each condition
cat results/acLDL/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz --bed results/acLDL/fastqtl/input/Ctrl.cqn_expression.txt.gz --cov results/acLDL/fastqtl/input/Ctrl.covariates.txt --W 500000 --permute '100 10000' --out results/acLDL/fastqtl/output/Ctrl_perm --execute True"
cat results/acLDL/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz --bed results/acLDL/fastqtl/input/AcLDL.cqn_expression.txt.gz --cov results/acLDL/fastqtl/input/AcLDL.covariates.txt --W 500000 --permute '100 10000' --out results/acLDL/fastqtl/output/AcLDL_perm --execute True"

#Merge chunks into single files
zcat results/acLDL/fastqtl/output/Ctrl_perm.chunk_*.txt.gz | bgzip > results/acLDL/fastqtl/output/Ctrl_permuted.txt.gz
zcat results/acLDL/fastqtl/output/AcLDL_perm.chunk_*.txt.gz | bgzip > results/acLDL/fastqtl/output/AcLDL_permuted.txt.gz

#Run on all SNPs 
cat results/acLDL/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz  --bed results/acLDL/fastqtl/input/Ctrl.cqn_expression.txt.gz --cov results/acLDL/fastqtl/input/Ctrl.covariates.txt --W 500000 --out results/acLDL/fastqtl/output/Ctrl_full --execute True"
cat results/acLDL/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz  --bed results/acLDL/fastqtl/input/AcLDL.cqn_expression.txt.gz --cov results/acLDL/fastqtl/input/AcLDL.covariates.txt --W 500000 --out results/acLDL/fastqtl/output/AcLDL_full --execute True"

#Merge chunks into a single file
zcat results/acLDL/fastqtl/output/Ctrl_full.chunk_*.txt.gz | bgzip > results/acLDL/fastqtl/output/Ctrl_pvalues.txt.gz
zcat results/acLDL/fastqtl/output/AcLDL_full.chunk_*.txt.gz | bgzip > results/acLDL/fastqtl/output/AcLDL_pvalues.txt.gz

#Remove the chunks
rm results/acLDL/fastqtl/output/Ctrl_perm.chunk_*
rm results/acLDL/fastqtl/output/AcLDL_perm.chunk_*

rm results/acLDL/fastqtl/output/Ctrl_full.chunk_*.txt.gz 
rm results/acLDL/fastqtl/output/AcLDL_full.chunk_*.txt.gz 

#Add SNP coordinates
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_07.vcf.gz --fastqtl results/acLDL/fastqtl/output/Ctrl_pvalues.txt.gz | bgzip > results/acLDL/fastqtl/output/Ctrl_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname run_fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/acLDL/fastqtl/input/fastqtl_genotypes.INFO_07.vcf.gz --fastqtl results/acLDL/fastqtl/output/AcLDL_pvalues.txt.gz | bgzip > results/acLDL/fastqtl/output/AcLDL_pvalues.coords.txt.gz"

#Sort files by SNP coordinates
#awk command is necessary to change field separator from space to tab
zcat results/acLDL/fastqtl/output/Ctrl_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/acLDL/fastqtl/output/Ctrl_pvalues.sorted.txt.gz &
zcat results/acLDL/fastqtl/output/AcLDL_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/acLDL/fastqtl/output/AcLDL_pvalues.sorted.txt.gz &


###### Map splicing QTLs using LeafCutter ######

#Construct junction files
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_sample_gt_map.txt | head -n 75 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamToJunctions --command "python ~/software/utils/bam/bamToJunctions.py --indir STAR/acLDL/ --outdir STAR/acLDL/ --leafCutterDir ~/software/leafcutter/ --insuffix .Aligned.sortedByCoord.RG.bam --execute True" 
echo "BEZI_24h_Ctrl" | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 2 --jobname bamToJunctions --command "python ~/software/utils/bam/bamToJunctions.py --indir STAR/acLDL/ --outdir STAR/acLDL/ --leafCutterDir ~/software/leafcutter/ --insuffix .Aligned.sortedByCoord.RG.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_sample_gt_map.txt | tail -n 75 | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 2 --jobname bamToJunctions --command "python ~/software/utils/bam/bamToJunctions.py --indir STAR/acLDL/ --outdir STAR/acLDL/ --leafCutterDir ~/software/leafcutter/ --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Create a list of junction files
ls --color=never STAR/acLDL/*/*.junc | cat > results/acLDL/leafcutter/junction_files.txt

#Cluster junctions
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o FarmOut/leafcutter_cluster.%J.jobout "python ~/software/leafcutter/clustering/leafcutter_cluster.py -j results/acLDL/leafcutter/junction_files.txt -r results/acLDL/leafcutter/ -m 50 -l 500000"

#Split leafcutter output into cluster counts and intron counts
python ~/software/utils/leafcutter/leafcutter_split_counts.py --leafcutter_out results/acLDL/leafcutter/leafcutter_perind.counts.gz --outprefix leafcutter

#Compress and index intron inclusion files
bgzip results/acLDL/leafcutter/fastqtl_input/Ctrl.norm_prop.txt && tabix -p bed results/acLDL/leafcutter/fastqtl_input/Ctrl.norm_prop.txt.gz
bgzip results/acLDL/leafcutter/fastqtl_input/AcLDL.norm_prop.txt && tabix -p bed results/acLDL/leafcutter/fastqtl_input/AcLDL.norm_prop.txt.gz

#Run FastQTL
cat results/acLDL/leafcutter/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz --bed results/acLDL/leafcutter/fastqtl_input/Ctrl.norm_prop.txt.gz --cov results/acLDL/leafcutter/fastqtl_input/Ctrl.covariates_prop.txt --W 100000 --permute '100 10000' --out results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_perm --execute True"
cat results/acLDL/leafcutter/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz --bed results/acLDL/leafcutter/fastqtl_input/AcLDL.norm_prop.txt.gz --cov results/acLDL/leafcutter/fastqtl_input/AcLDL.covariates_prop.txt --W 100000 --permute '100 10000' --out results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_perm --execute True"

#Merge chunks into single files
zcat results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_perm.chunk_*.txt.gz | bgzip > results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_permuted.txt.gz
zcat results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_perm.chunk_*.txt.gz | bgzip > results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_permuted.txt.gz

#Remove chunks
rm results/acLDL/leafcutter/fastqtl_output/*.chunk_*

#Get full p-values
cat results/acLDL/leafcutter/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz --bed results/acLDL/leafcutter/fastqtl_input/Ctrl.norm_prop.txt.gz --cov results/acLDL/leafcutter/fastqtl_input/Ctrl.covariates_prop.txt --W 100000 --out results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_full --execute True"
cat results/acLDL/leafcutter/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz --bed results/acLDL/leafcutter/fastqtl_input/AcLDL.norm_prop.txt.gz --cov results/acLDL/leafcutter/fastqtl_input/AcLDL.covariates_prop.txt --W 100000 --out results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_full --execute True"

#Merge chunks into single files
zcat results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_full.chunk_*.txt.gz | bgzip > results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_pvalues.txt.gz
zcat results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_full.chunk_*.txt.gz | bgzip > results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_pvalues.txt.gz

#Remove chunks
rm results/acLDL/leafcutter/fastqtl_output/*.chunk_*

#Add SNP coordinates
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz --fastqtl results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_pvalues.txt.gz | bgzip > results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz --fastqtl results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_pvalues.txt.gz | bgzip > results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_pvalues.coords.txt.gz"

#Sort files by SNP coordinates
#awk command is necessary to change field separator from space to tab
zcat results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/acLDL/leafcutter/fastqtl_output/Ctrl_100kb_pvalues.sorted.txt.gz &
zcat results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/acLDL/leafcutter/fastqtl_output/AcLDL_100kb_pvalues.sorted.txt.gz &

#Remove intermediate files
rm results/acLDL/leafcutter/fastqtl_output/*_pvalues.txt*




#Run FastQTL on the fold-change phenotype
bgzip results/acLDL/fastqtl/input_FC/FC.fold_change.txt && tabix -p bed results/acLDL/fastqtl/input_FC/FC.fold_change.txt.gz

#Run FastQTL
cat results/acLDL/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.INFO_07.vcf.gz --bed results/acLDL/fastqtl/input_FC/FC.fold_change.txt.gz --cov results/acLDL/fastqtl/input_FC/FC.covariates_PC3.txt --W 500000 --permute '100 10000' --out results/acLDL/fastqtl/output_FC/FC_perm --execute True"
zcat results/acLDL/fastqtl/output_FC/FC_perm.chunk_*.txt.gz | bgzip > results/acLDL/fastqtl/output_FC/FC_permuted.txt.gz
rm results/acLDL/fastqtl/output_FC/*.chunk_*






