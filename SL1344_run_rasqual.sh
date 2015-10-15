#Add read group to the BAM files to make them work with GATK
cat macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamAddRG --command "python ~/software/utils/bam/bamAddRG.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --outsuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Index bams
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Use ASEReadCounter to count allele-specific expression
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | head -n1 | python ~/software/utils/submitJobs.py --MEM 2000 --jobname bamCountASE --command "python ~/software/utils/bam/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.59_samples.snps_only.INFO_08.vcf.gz --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_sample_gt_map.txt | tail -n251 | python ~/software/utils/submitJobs.py --MEM 2000 --jobname bamCountASE --command "python ~/software/utils/bam/bamCountASE.py --indir STAR/SL1344 --outdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/snps_only/auotsomes.snps_only.sorted.uniq.no_duplicates.vcf --execute True"

#Merge ASE read counts per condition into one file:
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/counts/mergeASECounts.py --sample_list rasqual/input/SL1344_sg_map_A.txt --indir STAR/SL1344 --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condA.txt"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/counts/mergeASECounts.py --sample_list rasqual/input/SL1344_sg_map_B.txt --indir STAR/SL1344 --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condB.txt"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/counts/mergeASECounts.py --sample_list rasqual/input/SL1344_sg_map_C.txt --indir STAR/SL1344 --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condC.txt"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>5000] rusage[mem=5000]" -q normal -M 5000 -o FarmOut/mergeASEcounts.%J.jobout "python ~/software/utils/counts/mergeASECounts.py --sample_list rasqual/input/SL1344_sg_map_D.txt --indir STAR/SL1344 --suffix .ASEcounts > results/SL1344/ASE/SL1344_ASE_counts_condD.txt"

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

