#Lift over imputed genotypes from GRCh37 to GRCh38
cut -f1 genotypes/vcf_file_list.txt | tail -n 22 | python ~/software/utils/submitJobs.py --MEM 5000 --jobname liftOverVCF --command "python ~/software/utils/vcf/liftoverVcfGenotypes.py --chrMapFwd macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt --chrMapRev macrophage-gxe-study/data/liftOver_genotypes/Hg38ToGRCh38_chromosome_map.txt --liftOver macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain --reference ../../annotations/hg38/hg38.fa --vcfSuffix .vcf.gz --indir genotypes/raw/gtarray/imputed_vcf/20150128_858samples/ --outdir genotypes/GRCh38/ --execute True"

#Filter genotypes
echo "hipsci.chr21.gtarray.HumanCoreExome-12_v1_0.imputed_phased.858_samples.20150128.genotypes.GRCh38.sorted" | python ~/software/utils/submitJobs.py --MEM 5000 --jobname filterVCF --command "python ~/software/utils/vcf/filterVcf.py --sampleList genotypes/SL1344/SL1344_genotype_names.txt --MAF 0.05 --indir genotypes/GRCh38/imputed/ --outdir genotypes/SL1344/ --execute True --IMP2 0.7"

#Filter VCF files for samples, MAF and IMP2 values
cut -f1 genotypes/vcf_file_list.txt | head -n 22 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname filterVCF --command "python ~/software/utils/vcf/filterVcf.py  --sampleList macrophage-gxe-study/data/sample_lists/SL1344/SL1344_gt_list.txt --MAF 0.05 --indir genotypes/GRCh38/imputed/ --outdir genotypes/SL1344/imputed_filtered/ --execute True --IMP2 0.7 --vcfSuffix .GRCh38.sorted.vcf.gz"

#Filter VCF files again to only keep SNPs
cut -f1 genotypes/vcf_file_list.txt | head -n 22| python ~/software/utils/submitJobs.py --MEM 1000 --jobname filterVCFKeepSnps --command "python ~/software/utils/vcf/filterVcfKeepSnps.py  --indir genotypes/SL1344/imputed_filtered/ --outdir genotypes/SL1344/snps_only/ --execute True --vcfSuffix .filtered.GRCh38.sorted.vcf.gz"

#Concat chromosomes to count reads only once
# Note: this requires the indicudal vcf files to be indexed
bcftools concat -a chr1.snps_only.vcf.gz chr2.snps_only.vcf.gz chr3.snps_only.vcf.gz chr4.snps_only.vcf.gz chr5.snps_only.vcf.gz chr6.snps_only.vcf.gz chr7.snps_only.vcf.gz chr8.snps_only.vcf.gz chr9.snps_only.vcf.gz chr10.snps_only.vcf.gz chr11.snps_only.vcf.gz chr12.snps_only.vcf.gz chr13.snps_only.vcf.gz chr14.snps_only.vcf.gz chr15.snps_only.vcf.gz chr16.snps_only.vcf.gz chr17.snps_only.vcf.gz chr18.snps_only.vcf.gz chr19.snps_only.vcf.gz chr20.snps_only.vcf.gz chr21.snps_only.vcf.gz chr22.snps_only.vcf.gz > auotsomes.snps_only.vcf

#Sort genotypes (some chr1 snps went to chr9, GATK complains)
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o sortvcf.%J.jobout "~/software/vcflib/bin/vcfsort auotsomes.snps_only.vcf > auotsomes.snps_only.sorted.vcf"

#Run through vcfuniq, because liftover might have led to duplicated entries
~/software/vcflib/bin/vcfuniq auotsomes.snps_only.sorted.vcf > auotsomes.snps_only.sorted.uniq.vcf

#Sometimes there ALT allele of measured and genotyped variants do not agree. One of the variants has to be removed manually:
python ~/software/utils/vcf/vcfDetectDuplicateVariants.py --vcf auotsomes.snps_only.sorted.uniq.vcf > auotsomes.snps_only.sorted.uniq.no_duplicates.vcf

#Add read group to the BAM files to make them work with GATK
cat macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | head -n 132 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamAddRG --command "python ~/software/utils/bamAddRG.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --outsuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Index bams
cut -f1 genotypes/SL1344/SL1344_sample_genotype_map.txt | tail -n 108 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Use ASEReadCounter to count allele-specific expression
cut -f1 genotypes/SL1344/SL1344_sample_genotype_map.txt | head -n1 | python ~/software/utils/submitJobs.py --MEM 2000 --jobname bamCountASE --command "python ~/software/utils/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/snps_only/auotsomes.snps_only.sorted.uniq.no_duplicates.vcf --execute True"

#Merge ASE read counts per condition into one file:
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/mergeASECounts.py --sample_list genotypes/SL1344/SL1344_sg_map_A.txt --indir STAR/SL1344_ASE_counts --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condA.txt"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/mergeASECounts.py --sample_list genotypes/SL1344/SL1344_sg_map_B.txt --indir STAR/SL1344_ASE_counts --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condB.txt"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/mergeASECounts.py --sample_list genotypes/SL1344/SL1344_sg_map_C.txt --indir STAR/SL1344_ASE_counts --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condC.txt"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/mergeASECounts.py --sample_list genotypes/SL1344/SL1344_sg_map_D.txt --indir STAR/SL1344_ASE_counts --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condD.txt"

#Add ASE counts to the VCF file
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/vcfAddASE.%J.jobout "python ~/software/utils/vcf/vcfAddASE.py --ASEcounts results/SL1344/ASE/SL1344_ASE_counts_condA.txt --VCFfile genotypes/SL1344/snps_only/auotsomes.snps_only.sorted.uniq.no_duplicates.vcf | bgzip > results/SL1344/ASE/SL1344_ASE_counts_condA.vcf.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/vcfAddASE.%J.jobout "python ~/software/utils/vcf/vcfAddASE.py --ASEcounts results/SL1344/ASE/SL1344_ASE_counts_condB.txt --VCFfile genotypes/SL1344/snps_only/auotsomes.snps_only.sorted.uniq.no_duplicates.vcf | bgzip > results/SL1344/ASE/SL1344_ASE_counts_condB.vcf.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/vcfAddASE.%J.jobout "python ~/software/utils/vcf/vcfAddASE.py --ASEcounts results/SL1344/ASE/SL1344_ASE_counts_condC.txt --VCFfile genotypes/SL1344/snps_only/auotsomes.snps_only.sorted.uniq.no_duplicates.vcf | bgzip > results/SL1344/ASE/SL1344_ASE_counts_condC.vcf.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/vcfAddASE.%J.jobout "python ~/software/utils/vcf/vcfAddASE.py --ASEcounts results/SL1344/ASE/SL1344_ASE_counts_condD.txt --VCFfile genotypes/SL1344/snps_only/auotsomes.snps_only.sorted.uniq.no_duplicates.vcf | bgzip > results/SL1344/ASE/SL1344_ASE_counts_condD.vcf.gz"

#Extract exon start-end coordinates from the Txdb
bsub -G team170 -n1 -R "span[hosts=1] select[mem>12000] rusage[mem=12000]" -q normal -M 12000 -o FarmOut/exonCoords.%J.jobout "/software/R-3.1.2/bin/Rscript macrophage-gxe-study/SL1344/eQTL/convertTxdbToCoords.R"