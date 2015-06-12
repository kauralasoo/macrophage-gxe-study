#Fetch all file names from iRODS
python ~/software/utils/irodsGetSamplesInStudy.py --studyName "Genetics of gene expression in human macrophage response to Salmonella" |  cut -f1 -d "." | uniq > fastq/SL1344_samples.txt

#Match file names to sample names
python ~/software/utils/irodsFetchMeta.py --irodsList fastq/SL1344_samples.txt | sort -k1 > fastq/SL1344_names.txt 

#Download fastq files from iRODS
cut -f1 fastq/SL1344_samples.txt | head -n 100 |  python ~/software/utils/irodsToFastq.py --output fastq/SL1344/
cut -f1 fastq/SL1344_samples.txt | head -n 200 | tail -n 100 |  python ~/software/utils/irodsToFastq.py --output fastq/SL1344/
cut -f1 fastq/SL1344_samples.txt | tail -n 328 |  python ~/software/utils/irodsToFastq.py --output fastq/SL1344/

#Convert remaining crams to fastq
cut -f1 fastq/remaining.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/cramToBam.py --inputDir fastq/SL1344/ --outputDir fastq/SL1344/"

#Merge split bams into joint ones and rename
cat fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/ --outdir fastq/SL1344/ --suffix .1.fastq.gz"
cat fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/ --outdir fastq/SL1344/ --suffix .2.fastq.gz"

#Align reads to the transcriptome using STAR
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR/SL1344/ --fastqDir fastq/SL1344/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"

cut -f1 fastq/SL1344_names.txt | tail -n 16 | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR/SL1344/ --fastqDir fastq/SL1344/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"

#Rerun failed samples
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR/SL1344_rerun/ --fastqDir fastq/SL1344/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"


#Count the number of reads overlapping gene annotations
cut -f1 fastq/SL1344_names.txt | tail -n 40 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.78.gtf --strand 2 --execute True"

cut -f1 fastq/SL1344_names.txt | head -n 92 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.78.gtf --strand 2 --execute True"

cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname mycoplasmaTest --command "python ~/software/utils/mycoplasmaTest.py --inputDir fastq/SL1344/ --outdir STAR/SL1344/ --bwaIndex ../../annotations/Mycoplasma/bwa_index/Mycoplasma_genomes.fa --execute True"

cut -f1 fastq/SL1344_names.txt | head -n 1 | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR1 --fastqDir fastq/SL1344/ --genomeDir ../../annotations/GRCh38/STAR_index/ --runThreadN 8"

#Count reads overlaping GENCODE basic annotations
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"

#Count reads overlapping full Ensembl 79 annotations
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gtf --strand 2 --countsSuffix .counts.txt --execute True"
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gtf --strand 2 --countsSuffix .counts.txt --execute True"

#Convert bedgraph files to bigwig
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"

cut -f1 fastq/SL1344_names.txt | tail -n 16 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 fastq/SL1344_names.txt | tail -n 16 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"

#Index bam files
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --execute True"


#Use verifyBamID to check concordance with the vcf file
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/selected_genotypes.GRCh38.sorted.vcf.gz --execute True" 
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/selected_genotypes.GRCh38.sorted.vcf.gz --execute True" 

#Remove SNPs with multiple IDs
grep -v ";" eqtlbma/snp_coords.bed > eqtlbma/snp_coords_filtered.bed
cat eqtlbma/snp_coords_filtered.bed | sort -k1,1V -k2,2g | bgzip > eqtlbma/snp_coords_filtered.bed.gz
tabix -p bed eqtlbma/snp_coords_filtered.bed.gz 


#Run eqtlbma on a small number of genes
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o FarmOut/eqtlbma_bf.%J.jobout "eqtlbma_bf --geno eqtlbma/list_genotypes.txt --scoord eqtlbma/snp_coords.bed.gz --exp eqtlbma/list_explevels.txt --gcoord eqtlbma/gene_coords.subset.sort.bed --anchor TSS --cis 500000 --out eqtlbma/output/subset_test --analys join --gridL eqtlbma/grid_phi2_oma2_general.txt  --gridS eqtlbma/grid_phi2_oma2_with-configs.txt --bfs all --error mvlr -v 3"

bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o FarmOut/eqtlbma_hm.%J.jobout  "eqtlbma_hm --data 'eqtlbma/output/subset_test_l10abfs_raw.txt.gz' --nsubgrp 4 --dim 15 --ngrid 10 --out eqtlbma/output/subset_test_eqtlbma_hm_txt.gz"

#Extract config weights
zcat eqtlbma/output/subset_test_eqtlbma_hm_txt.gz | grep "#config" | awk '{split($1,a,"."); print a[2]"\t"$2}' > eqtlbma/output/subset_test_config_weights.txt

#Extract grid weigths
zcat eqtlbma/output/subset_test_eqtlbma_hm_txt.gz | grep "#grid" | cut -f2 > eqtlbma/output/subset_test_grid_weights.txt

#Lift over imputed genotypes from GRCh37 to GRCh38
cut -f1 genotypes/vcf_file_list.txt | tail -n 22 | python ~/software/utils/submitJobs.py --MEM 5000 --jobname liftOverVCF --command "python ~/software/utils/vcf/liftoverVcfGenotypes.py --chrMapFwd macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt --chrMapRev macrophage-gxe-study/data/liftOver_genotypes/Hg38ToGRCh38_chromosome_map.txt --liftOver macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain --reference ../../annotations/hg38/hg38.fa --vcfSuffix .vcf.gz --indir genotypes/raw/gtarray/imputed_vcf/20150128_858samples/ --outdir genotypes/GRCh38/ --execute True"

#LiftOver raw genotypes from GRCh37 to GRCh38
echo "hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes" | python ~/software/utils/submitJobs.py --MEM 5000 --jobname liftOverVCF --command "python ~/software/utils/vcf/liftoverVcfGenotypes.py --chrMapFwd macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt --chrMapRev macrophage-gxe-study/data/liftOver_genotypes/Hg38ToGRCh38_chromosome_map.txt --liftOver macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain --reference ../../annotations/hg38/hg38.fa --vcfSuffix .vcf.gz --indir genotypes/raw/gtarray/vcf/20141111_858samples/ --outdir genotypes/GRCh38/ --execute True"


cut -f1 genotypes/vcf_file_list.txt | grep chr3 | python ~/software/utils/submitJobs.py --MEM 5000 --jobname liftOverVCF --command "python ~/software/utils/vcf/liftoverVcfGenotypes.py --chrMapFwd macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt --chrMapRev macrophage-gxe-study/data/liftOver_genotypes/Hg38ToGRCh38_chromosome_map.txt --liftOver macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain --reference ../../annotations/hg38/hg38.fa --vcfSuffix .vcf.gz --indir genotypes/raw/gtarray/imputed_vcf/20150128_858samples/ --outdir genotypes/GRCh38/ --execute True"

cut -f1 genotypes/vcf_file_list.txt | grep chr6 | python ~/software/utils/submitJobs.py --MEM 5000 --jobname liftOverVCF --command "python ~/software/utils/vcf/liftoverVcfGenotypes.py --chrMapFwd macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt --chrMapRev macrophage-gxe-study/data/liftOver_genotypes/Hg38ToGRCh38_chromosome_map.txt --liftOver macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain --reference ../../annotations/hg38/hg38.fa --vcfSuffix .vcf.gz --indir genotypes/raw/gtarray/imputed_vcf/20150128_858samples/ --outdir genotypes/GRCh38/ --execute True"

cut -f1 genotypes/vcf_file_list.txt | grep chr7 | python ~/software/utils/submitJobs.py --MEM 5000 --jobname liftOverVCF --command "python ~/software/utils/vcf/liftoverVcfGenotypes.py --chrMapFwd macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt --chrMapRev macrophage-gxe-study/data/liftOver_genotypes/Hg38ToGRCh38_chromosome_map.txt --liftOver macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain --reference ../../annotations/hg38/hg38.fa --vcfSuffix .vcf.gz --indir genotypes/raw/gtarray/imputed_vcf/20150128_858samples/ --outdir genotypes/GRCh38/ --execute True"

#Filter genotypes
echo "hipsci.chr21.gtarray.HumanCoreExome-12_v1_0.imputed_phased.858_samples.20150128.genotypes.GRCh38.sorted" | python ~/software/utils/submitJobs.py --MEM 5000 --jobname filterVCF --command "python ~/software/utils/vcf/filterVcf.py --sampleList genotypes/SL1344/SL1344_genotype_names.txt --MAF 0.05 --indir genotypes/GRCh38/imputed/ --outdir genotypes/SL1344/ --execute True --IMP2 0.7"


