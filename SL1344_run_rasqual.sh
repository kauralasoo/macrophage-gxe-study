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
python ~/software/utils/vcf/vcfDetectDuplicateVariants.py --vcf auotsomes.snps_only.sorted.uniq.vcf --duplicates duplicates.txt > auotsomes.snps_only.sorted.uniq.no_duplicates.vcf

#Add read group to the BAM files to make them work with GATK
cat macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | head -n 132 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamAddRG --command "python ~/software/utils/bamAddRG.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --outsuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Index bams
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Use ASEReadCounter to count allele-specific expression
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | head -n1 | python ~/software/utils/submitJobs.py --MEM 2000 --jobname bamCountASE --command "python ~/software/utils/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/snps_only/auotsomes.snps_only.sorted.uniq.no_duplicates.vcf --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | tail -n251 | python ~/software/utils/submitJobs.py --MEM 2000 --jobname bamCountASE --command "python ~/software/utils/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/snps_only/auotsomes.snps_only.sorted.uniq.no_duplicates.vcf --execute True"

#Merge ASE read counts per condition into one file:
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/mergeASECounts.py --sample_list rasqual/input/SL1344_sg_map_A.txt --indir STAR/SL1344 --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condA.txt"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/mergeASECounts.py --sample_list rasqual/input/SL1344_sg_map_B.txt --indir STAR/SL1344 --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condB.txt"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/mergeASECounts.py --sample_list rasqual/input/SL1344_sg_map_C.txt --indir STAR/SL1344 --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condC.txt"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/mergeASECounts.py --sample_list rasqual/input/SL1344_sg_map_D.txt --indir STAR/SL1344 --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condD.txt"

#Construct genotype list
cut -f1 SL1344_sg_map_A.txt > genotype_list.txt 

#Extract only relevant genotypes from the vcf file (in the correct order)
bcftools view -S rasqual/input/genotype_list.txt genotypes/SL1344/snps_only/auotsomes.snps_only.sorted.uniq.no_duplicates.vcf > genotypes/SL1344/autosomes.snps_only.no_replicates.vcf
echo "hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname filterVCF --command "python ~/software/utils/vcf/filterVcf.py  --sampleList rasqual/input/genotype_list.txt --MAF 0.05 --indir genotypes/GRCh38/genotyped/ --outdir genotypes/SL1344/ --execute True  --vcfSuffix .vcf.gz"

#Add ASE counts to the VCF file
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/vcfAddASE.%J.jobout "python ~/software/utils/vcf/vcfAddASE.py --ASEcounts results/SL1344/ASE/SL1344_ASE_counts_condA.txt --VCFfile genotypes/SL1344/autosomes.snps_only.no_replicates.vcf | bgzip > rasqual/input/ASE/SL1344_ASE_counts_condA.vcf.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/vcfAddASE.%J.jobout "python ~/software/utils/vcf/vcfAddASE.py --ASEcounts results/SL1344/ASE/SL1344_ASE_counts_condB.txt --VCFfile genotypes/SL1344/autosomes.snps_only.no_replicates.vcf | bgzip > rasqual/input/ASE/SL1344_ASE_counts_condB.vcf.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/vcfAddASE.%J.jobout "python ~/software/utils/vcf/vcfAddASE.py --ASEcounts results/SL1344/ASE/SL1344_ASE_counts_condC.txt --VCFfile genotypes/SL1344/autosomes.snps_only.no_replicates.vcf | bgzip > rasqual/input/ASE/SL1344_ASE_counts_condC.vcf.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/vcfAddASE.%J.jobout "python ~/software/utils/vcf/vcfAddASE.py --ASEcounts results/SL1344/ASE/SL1344_ASE_counts_condD.txt --VCFfile genotypes/SL1344/autosomes.snps_only.no_replicates.vcf | bgzip > rasqual/input/ASE/SL1344_ASE_counts_condD.vcf.gz"

#Index the VCF files using tabix
tabix -p vcf rasqual/input/ASE/SL1344_ASE_counts_condA.vcf.gz
tabix -p vcf rasqual/input/ASE/SL1344_ASE_counts_condB.vcf.gz
tabix -p vcf rasqual/input/ASE/SL1344_ASE_counts_condC.vcf.gz
tabix -p vcf rasqual/input/ASE/SL1344_ASE_counts_condD.vcf.gz

#Extract exon start-end coordinates from the Txdb
bsub -G team170 -n1 -R "span[hosts=1] select[mem>12000] rusage[mem=12000]" -q normal -M 12000 -o FarmOut/exonCoords.%J.jobout "/software/R-3.1.2/bin/Rscript macrophage-gxe-study/SL1344/eQTL/convertTxdbToCoords.R"

#Extract SNP coords from a vcf file
grep -v "#" genotypes/SL1344/autosomes.snps_only.no_replicates.vcf | cut -f 1,2,3 > genotypes/SL1344/snp_coords.txt

#Count SnpPerGene
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/countSnpsPerGene.%J.jobout "/software/R-3.1.2/bin/Rscript ~/software/utils/vcf/countsSnpsPerGene.R -s genotypes/SL1344/snp_coords.txt -e annotations/Homo_sapiens.GRCh38.79.gene_exon_start_end.filtered_genes.txt -c genotypes/SL1344/snps_per_gene.txt"


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

