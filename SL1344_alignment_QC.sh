#Fetch all file names from iRODS
python ~/software/utils/irodsGetSamplesInStudy.py --studyName "Genetics of gene expression in human macrophage response to Salmonella" |  cut -f1 -d "." | uniq > fastq/SL1344_samples3.txt

#Match file names to sample names
python ~/software/utils/irodsFetchMeta.py --irodsList fastq/SL1344_samples.txt | sort -k1 > fastq/SL1344_names.txt 
python ~/software/utils/irodsFetchMeta.py --irodsList fastq/SL1344_samples_2.txt | sort -k1 > fastq/SL1344_names_2.txt 

#Fetch lanelets in cram format from irods
cut -f1 fastq/acLDL_samples.txt | python ~/software/utils/fetch-irods.py --dir fastq/acLDL/cram/ --suffix .cram

#Convert remaining crams to fastq
cut -f1 fastq/remaining.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/cramToBam.py --inputDir fastq/SL1344/ --outputDir fastq/SL1344/"

#Merge split bams into joint ones and rename
cat fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/ --outdir fastq/SL1344/ --suffix .1.fastq.gz"
cat fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/ --outdir fastq/SL1344/ --suffix .2.fastq.gz"

#Align reads to the transcriptome using STAR
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR/SL1344/ --fastqDir fastq/SL1344/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"

#Count the number of reads aligning to the mycoplasma genome
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname mycoplasmaTest --command "python ~/software/utils/mycoplasmaTest.py --inputDir fastq/SL1344/ --outdir STAR/SL1344/ --bwaIndex ../../annotations/Mycoplasma/bwa_index/Mycoplasma_genomes.fa --execute True"

#Count reads overlaping GENCODE basic annotations
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"

#Count reads overlapping full Ensembl 79 annotations
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gtf --strand 2 --countsSuffix .counts.txt --execute True"

#Convert bedgraph files to bigwig
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

#Have to run this again!!!!
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

#Index bam files
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --execute True"

#LiftOver raw genotypes from GRCh37 to GRCh38
echo "hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes" | python ~/software/utils/submitJobs.py --MEM 5000 --jobname liftOverVCF --command "python ~/software/utils/vcf/liftoverVcfGenotypes.py --chrMapFwd macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt --chrMapRev macrophage-gxe-study/data/liftOver_genotypes/Hg38ToGRCh38_chromosome_map.txt --liftOver macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain --reference ../../annotations/hg38/hg38.fa --vcfSuffix .vcf.gz --indir genotypes/raw/gtarray/vcf/20141111_858samples/ --outdir genotypes/GRCh38/ --execute True"

#Extract samples from the large VCF file
echo "hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname filterVCF --command "python ~/software/utils/vcf/filterVcf.py  --sampleList genotypes/acLDL/acLDL_genotype_list.txt --MAF 0.05 --indir genotypes/GRCh38/genotyped/ --outdir genotypes/acLDL/ --execute True  --vcfSuffix .vcf.gz"

#Use verifyBamID to check concordance with the vcf file
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/selected_genotypes.GRCh38.sorted.vcf.gz --execute True" 
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/selected_genotypes.GRCh38.sorted.vcf.gz --execute True" 



