#Add read group to the BAM files to make them work with GATK
cat macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamAddRG --command "python ~/software/utils/bam/bamAddRG.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --outsuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Index bams
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Index the vcf file
tabix -p vcf genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz
bcftools index genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz

#Use ASEReadCounter to count allele-specific expression
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | python ~/software/utils/submitJobs.py --MEM 6000 --jobname bamCountASE --command "python ~/software/utils/rasqual/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True --Xmx 4g"

cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | tail -n251 | python ~/software/utils/submitJobs.py --MEM 2000 --jobname bamCountASE --command "python ~/software/utils/rasqual/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | head -n 45 | tail -n 44 | python ~/software/utils/submitJobs.py --MEM 2000 --jobname bamCountASE --command "python ~/software/utils/rasqual/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True"
echo "pelm_A" | python ~/software/utils/submitJobs.py --MEM 6000 --jobname bamCountASE --command "python ~/software/utils/rasqual/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True --Xmx 4g"

echo "pelm_A" | python ~/software/utils/rasqual/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute False --Xmx 4g

echo "jorr_A" | python ~/software/utils/submitJobs.py --MEM 10000 --jobname bamCountASE --command "python ~/software/utils/rasqual/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True --Xmx 8g"

#Construct a sample-sample map for meregeASECounts.py script
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | awk -v OFS='\t' '{print $1, $1}' > results/SL1344/rasqual/input/sample_sample_map.txt

#Merge all allele-specific counts into one matrix
echo "mergeASECounts" | python ~/software/utils/submitJobs.py --MEM 18000 --jobname mergeASECounts --command "python ~/software/utils/rasqual/mergeASECounts.py --sample_list results/SL1344/rasqual/input/sample_sample_map.txt --indir STAR/SL1344 --suffix .ASEcounts > results/SL1344/combined_ASE_counts.txt"

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
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 32000 --jobname vcfAddASE --queue hugemem --command "python ~/software/utils/rasqual/vcfAddASE.py --ASEcounts results/SL1344/combined_ASE_counts.txt --ASESampleGenotypeMap results/SL1344/rasqual/input/naive.sg_map.txt --VCFfile results/SL1344/rasqual/input/naive.vcf | bgzip > results/SL1344/rasqual/input/naive.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 32000 --jobname vcfAddASE --queue hugemem --command "python ~/software/utils/rasqual/vcfAddASE.py --ASEcounts results/SL1344/combined_ASE_counts.txt --ASESampleGenotypeMap results/SL1344/rasqual/input/IFNg.sg_map.txt --VCFfile results/SL1344/rasqual/input/IFNg.vcf | bgzip > results/SL1344/rasqual/input/IFNg.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 32000 --jobname vcfAddASE --queue hugemem --command "python ~/software/utils/rasqual/vcfAddASE.py --ASEcounts results/SL1344/combined_ASE_counts.txt --ASESampleGenotypeMap results/SL1344/rasqual/input/SL1344.sg_map.txt --VCFfile results/SL1344/rasqual/input/SL1344.vcf | bgzip > results/SL1344/rasqual/input/SL1344.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 32000 --jobname vcfAddASE --queue hugemem --command "python ~/software/utils/rasqual/vcfAddASE.py --ASEcounts results/SL1344/combined_ASE_counts.txt --ASESampleGenotypeMap results/SL1344/rasqual/input/IFNg_SL1344.sg_map.txt --VCFfile results/SL1344/rasqual/input/IFNg_SL1344.vcf | bgzip > results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz"

#Index VCF files using tabix
tabix -p vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz &
tabix -p vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz &
tabix -p vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz &
tabix -p vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz &

#Extract exon start-end coordinates from the Txdb
bsub -G team170 -n1 -R "span[hosts=1] select[mem>12000] rusage[mem=12000]" -q normal -M 12000 -o FarmOut/exonCoords.%J.jobout "/software/R-3.1.2/bin/Rscript macrophage-gxe-study/SL1344/eQTL/convertTxdbToCoords.R"

#RASQUAL (2 PEER covariates + sex, library sizes + GC), 500kb
cat results/SL1344/rasqual/input/gene_batches.txt | head -n 201 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_naive_500kb --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/naive_500kb --covariates results/SL1344/rasqual/input/naive.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force' --execute True"

cat results/SL1344/rasqual/input/gene_batches.txt | head -n 1 | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_IFNg_500kb --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/IFNg.expression.bin  --offsets results/SL1344/rasqual/input/IFNg.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/IFNg_500kb --covariates results/SL1344/rasqual/input/IFNg.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force' --execute True"
cat results/SL1344/rasqual/input/gene_batches.txt | tail -n 1200 | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_IFNg_500kb --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/IFNg.expression.bin  --offsets results/SL1344/rasqual/input/IFNg.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/IFNg_500kb --covariates results/SL1344/rasqual/input/IFNg.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force' --execute True"

cat results/SL1344/rasqual/input/gene_batches.txt | head -n 1 | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_SL1344_500kb --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/SL1344.expression.bin  --offsets results/SL1344/rasqual/input/SL1344.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/SL1344_500kb --covariates results/SL1344/rasqual/input/SL1344.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force' --execute True"
cat results/SL1344/rasqual/input/gene_batches.txt | tail -n 1200 | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_SL1344_500kb --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/SL1344.expression.bin  --offsets results/SL1344/rasqual/input/SL1344.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/SL1344_500kb --covariates results/SL1344/rasqual/input/SL1344.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force' --execute True"


cat results/SL1344/rasqual/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_IFNg_SL1344_500kb --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/IFNg_SL1344.expression.bin  --offsets results/SL1344/rasqual/input/IFNg_SL1344.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_500kb.txt --outprefix results/SL1344/rasqual/output/IFNg_SL1344_500kb --covariates results/SL1344/rasqual/input/IFNg_SL1344.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force' --execute True"


#Run RASQUAL on a subset of genes
echo "ENSG00000170458,ENSG00000166750,ENSG00000104689,ENSG00000109861" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command  "python ~/software/utils/runRasqual.py --readCounts rasqual/input/counts/cond_A_counts.bin --offsets rasqual/input/counts/cond_A_factors.bin --n 59 --geneids rasqual/input/counts/gene_id_list.txt --vcf rasqual/input/ASE/SL1344_ASE_counts_condA.vcf.gz --geneMetadata rasqual/input/snp_counts.txt --featureCoords annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.txt > rasqual.out"

echo "ENSG00000170458" | python ~/software/utils/runRasqual.py --readCounts rasqual/input/counts/cond_A_counts.bin --offsets rasqual/input/counts/cond_A_factors.bin --n 59 --geneids rasqual/input/counts/gene_id_list.txt --vcf rasqual/input/ASE/SL1344_ASE_counts_condA.vcf.gz --geneMetadata genotypes/SL1344/snps_per_gene.txt --execute True > rasqual.out3

#VCAM1
echo "ENSG00000162692" | python ~/software/utils/runRasqual.py --readCounts rasqual/input/counts/cond_A_counts.bin --offsets rasqual/input/counts/cond_A_factors.bin --n 59 --geneids rasqual/input/counts/gene_id_list.txt --vcf rasqual/input/ASE/SL1344_ASE_counts_condA.vcf.gz --geneMetadata genotypes/SL1344/snps_per_gene.txt --execute True > VCAM1.cond_A.rasqual
echo "ENSG00000162692" | python ~/software/utils/runRasqual.py --readCounts rasqual/input/counts/cond_B_counts.bin --offsets rasqual/input/counts/cond_B_factors.bin --n 59 --geneids rasqual/input/counts/gene_id_list.txt --vcf rasqual/input/ASE/SL1344_ASE_counts_condB.vcf.gz --geneMetadata genotypes/SL1344/snps_per_gene.txt --execute True > VCAM1.cond_B.rasqual

python ~/software/utils/rasqualToIGV.py --rasqualOut VCAM1.cond_A.rasqual --geneid ENSG00000162692 > VCAM1.cond_A.gwas
python ~/software/utils/rasqualToIGV.py --rasqualOut VCAM1.cond_B.rasqual --geneid ENSG00000162692 > VCAM1.cond_B.gwas


#CTSC
echo "ENSG00000109861" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts rasqual/input/counts/cond_A_counts.bin --offsets rasqual/input/counts/cond_A_factors.bin --n 59 --geneids rasqual/input/counts/gene_id_list.txt --vcf rasqual/input/ASE/SL1344_ASE_counts_condA.vcf.gz --geneMetadata genotypes/SL1344/snps_per_gene.txt --execute True > CTSC.cond_A.rasqual"
echo "ENSG00000109861" |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts rasqual/input/counts/cond_B_counts.bin --offsets rasqual/input/counts/cond_B_factors.bin --n 59 --geneids rasqual/input/counts/gene_id_list.txt --vcf rasqual/input/ASE/SL1344_ASE_counts_condB.vcf.gz --geneMetadata genotypes/SL1344/snps_per_gene.txt --execute True > CTSC.cond_B.rasqual"

python ~/software/utils/rasqualToIGV.py --rasqualOut CTSC.cond_A.rasqual --geneid ENSG00000109861 > CTSC.cond_A.gwas
python ~/software/utils/rasqualToIGV.py --rasqualOut CTSC.cond_B.rasqual --geneid ENSG00000109861 > CTSC.cond_B.gwas



