#Fetch all file names from iRODS
python ~/software/utils/irods/irodsGetSamplesInStudy.py --studyName "Genetics of gene expression in human macrophage response to Salmonella" |  cut -f1 -d "." | uniq > ATAC_samples_5.txt

#Match irods file names to sample names
python ~/software/utils/irods/irodsFetchMeta.py --irodsList data/SL1344/ATAC_Salmonella_samples.txt | sort -k1 > data/SL1344/ATAC_Salmonella_names.txt 
python ~/software/utils/irods/irodsFetchMeta.py --irodsList data/SL1344/ATAC_Salmonella_samples_2.txt | sort -k1 > data/SL1344/ATAC_Salmonella_names_2.txt 
python ~/software/utils/irods/irodsFetchMeta.py --irodsList macrophage-chromatin/data/SL1344/ATAC_Salmonella_samples_3.txt | sort -k1 > macrophage-chromatin/data/SL1344/ATAC_Salmonella_sample_names_3.txt

#Fetch CRAM files from irods as well
cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_samples.txt | python ~/software/utils/irods/fetch-irods.py --dir data/ATAC_Salmonella/cram/ --suffix .cram
cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_samples_2.txt | python ~/software/utils/irods/fetch-irods.py --dir data/ATAC_Salmonella/cram/ --suffix .cram
cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_samples_3.txt | python ~/software/utils/irods/fetch-irods.py --dir data/ATAC_Salmonella/cram/ --suffix .cram

#Convert CRAMs to FASTQ
cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_samples.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/bam/cramToFastq.py --inputDir data/ATAC_Salmonella/cram/ --outputDir data/ATAC_Salmonella/fastq/"
cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_samples_2.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/bam/cramToFastq.py --inputDir data/ATAC_Salmonella/cram/ --outputDir data/ATAC_Salmonella/fastq/"
cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_samples_3.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/bam/cramToFastq.py --inputDir data/ATAC_Salmonella/cram/ --outputDir data/ATAC_Salmonella/fastq/"

#Merge split fastq into joint ones and rename
cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_fastq --command "python ~/software/utils/fastq/merge-fastq.py --indir data/ATAC_Salmonella/fastq/ --outdir data/ATAC_Salmonella/fastq/ --suffix .1.fastq.gz"
cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_fastq --command "python ~/software/utils/fastq/merge-fastq.py --indir data/ATAC_Salmonella/fastq/ --outdir data/ATAC_Salmonella/fastq/ --suffix .2.fastq.gz"

cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_fastq --command "python ~/software/utils/fastq/merge-fastq.py --indir data/ATAC_Salmonella/fastq/ --outdir data/ATAC_Salmonella/fastq/ --suffix .1.fastq.gz"
cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_fastq --command "python ~/software/utils/fastq/merge-fastq.py --indir data/ATAC_Salmonella/fastq/ --outdir data/ATAC_Salmonella/fastq/ --suffix .2.fastq.gz"

cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_fastq --command "python ~/software/utils//fastq/merge-fastq.py --indir data/ATAC_Salmonella/fastq/ --outdir data/ATAC_Salmonella/fastq/ --suffix .1.fastq.gz"
cat  macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_fastq --command "python ~/software/utils/fastq/merge-fastq.py --indir data/ATAC_Salmonella/fastq/ --outdir data/ATAC_Salmonella/fastq/ --suffix .2.fastq.gz"

### ALIGNMENT AND FILTERING ###
#Extract sample indexes from metadata
cut -f1,4,6 -d, --output-delimiter=$'\t' macrophage-chromatin/data/SL1344/ATAC_sample_metadata.csv > macrophage-chromatin/data/SL1344/ATAC_sample_indexes.txt

#Trim adapters using skewer
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | head -n 40 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname fastq_trim --command "python ~/software/utils/fastq/ATAC_trim_adapters.py --fastqDirIn data/ATAC_Salmonella/fastq_untrimmed --fastqDirOut data/ATAC_Salmonella/fastq_trimmed/ --atacPrimers macrophage-chromatin/data/SL1344/ATAC_index_primers.txt --sampleIndexMap macrophage-chromatin/data/SL1344/ATAC_sample_indexes.txt"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | tail -n 42 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname fastq_trim --command "python ~/software/utils/fastq/ATAC_trim_adapters.py --fastqDirIn data/ATAC_Salmonella/fastq_untrimmed --fastqDirOut data/ATAC_Salmonella/fastq_trimmed/ --atacPrimers macrophage-chromatin/data/SL1344/ATAC_index_primers.txt --sampleIndexMap macrophage-chromatin/data/SL1344/ATAC_sample_indexes.txt"

cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | grep "babk" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname fastq_trim --command "python ~/software/utils/fastq/ATAC_trim_adapters.py --fastqDirIn data/ATAC_Salmonella/fastq --fastqDirOut data/ATAC_Salmonella/fastq_trimmed/ --atacPrimers macrophage-chromatin/data/SL1344/ATAC_index_primers.txt --sampleIndexMap macrophage-chromatin/data/SL1344/ATAC_sample_indexes.txt"

#Use bwa to align the reads to the genome
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 20000 --ncores 8 --jobname bwa_mem --command  "python ~/software/utils/align/bwaAlignPE.py --outputDir processed/SL1344/ --fastqDir data/ATAC_Salmonella/fastq_trimmed/ --nCores 8 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38 --fastqSuffix .trimmed.fastq.gz"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | grep "babk" | python ~/software/utils/submitJobs.py --MEM 20000 --ncores 6 --jobname bwa_mem --command  "python ~/software/utils/align/bwaAlignPE.py --outputDir processed/SL1344/ --fastqDir data/ATAC_Salmonella/fastq_trimmed/ --nCores 6 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38 --fastqSuffix .trimmed.fastq.gz"

#Sort bams by coordinate 
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | head -n 20 | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/SL1344/ --outdir processed/SL1344/"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | head -n 50 | tail -n 30 | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/SL1344/ --outdir processed/SL1344/"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | tail -n 32 | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/SL1344/ --outdir processed/SL1344/"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | head -n 34 | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/SL1344/ --outdir processed/SL1344/"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | head -n 58 | tail -n 14 | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/SL1344/ --outdir processed/SL1344/"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | tail -n 10 | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/SL1344/ --outdir processed/SL1344/"


#Index bams
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/SL1344/ --insuffix .sortedByCoord.bam --execute True"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/SL1344/ --insuffix .sortedByCoord.bam --execute True"

#Keep only properly paired reads on primary chromosomes (remove MT and contigs)
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 1 --jobname remove_contigs --command "python ~/software/utils/bam/bamRemoveContigsMT.py --indir processed/SL1344/ --outdir processed/SL1344/"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 1 --jobname remove_contigs --command "python ~/software/utils/bam/bamRemoveContigsMT.py --indir processed/SL1344/ --outdir processed/SL1344/"

#Remove BWA header because it clashes with the Picard MarkDuplicates
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 1 --jobname remove_header --command "python ~/software/utils/bam/bamRemoveBwaHeader.py --indir processed/SL1344/ --outdir processed/SL1344/"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 1 --jobname remove_header --command "python ~/software/utils/bam/bamRemoveBwaHeader.py --indir processed/SL1344/ --outdir processed/SL1344/"

#Remove duplicates
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 1 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir processed/SL1344/ --outdir processed/SL1344/ --execute True --insuffix .reheadered.bam"
echo 'cicb_A_ATAC' | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 1 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir processed/SL1344/ --outdir processed/SL1344/ --execute True --insuffix .reheadered.bam"
echo 'oefg_B_ATAC' | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 4 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir processed/SL1344/ --outdir processed/SL1344/ --execute True --insuffix .reheadered.bam"
echo 'sukz_B_ATAC' | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 4 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir processed/SL1344/ --outdir processed/SL1344/ --execute True --insuffix .reheadered.bam"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 4 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir processed/SL1344/ --outdir processed/SL1344/ --execute True --insuffix .reheadered.bam"


#### QC and cleanup ####
#Count number of reads for each chromosome
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname count_chrs --command "python ~/software/utils/bam/countReadsPerChr.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .sortedByCoord.bam"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname count_chrs --command "python ~/software/utils/bam/countReadsPerChr.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .sortedByCoord.bam"

#Remove intermediate files
rm processed/SL1344/*/*.reheadered.bam
rm processed/SL1344/*/*.new_header.txt
rm processed/SL1344/*/*.sortedByCoord.bam*
rm processed/SL1344/*/*.filtered.bam
rm processed/SL1344/*/*.sortedByName.bam
ls processed/SL1344/*/*ATAC.bam

#Index bams
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/SL1344/ --insuffix .no_duplicates.bam --execute True"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/SL1344/ --insuffix .no_duplicates.bam --execute True"

#Extract genotypes from the chipped vcf file
#Only keep the relevant samples from the vcf file
bcftools view -S vcf_sample_names.txt ../macrophage-gxe-study/genotypes/GRCh38/genotyped/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted.vcf.gz | bcftools filter -O z -i 'MAF[0] >= 0.05' - > data/ATAC_Salmonella/genotypes/array_genotypes.86_samples.vcf.gz

#Run verifyBamID
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 1500 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir processed/SL1344/ --insuffix .no_duplicates.bam --vcf data/ATAC_Salmonella/genotypes/array_genotypes.27_samples.vcf.gz --execute True" 

cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1500 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir processed/SL1344/ --insuffix .no_duplicates.bam --vcf data/ATAC_Salmonella/genotypes/array_genotypes.86_samples.vcf.gz --execute True" 

#Run verifyBamID on Natsuhikos samples
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | grep "_D_" | python ~/software/utils/submitJobs.py --MEM 1500 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir ~/group-scratch/MacrophageATAC/Bams/ --insuffix .bam --vcf data/ATAC_Salmonella/genotypes/array_genotypes.86_samples.vcf.gz --execute True" 
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | grep "_A_" | python ~/software/utils/submitJobs.py --MEM 1500 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir ~/group-scratch/MacrophageATAC/Bams/ --insuffix .bam --vcf data/ATAC_Salmonella/genotypes/array_genotypes.86_samples.vcf.gz --execute True" 

#### Create fragment coverage files ####

#Sort bams by name
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | head -n42 | python ~/software/utils/submitJobs.py --MEM 9000 --ncores 1 --jobname bamsort_name --command "python ~/software/utils/bam/bamSortName.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .no_duplicates.bam"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | tail -n40 | python ~/software/utils/submitJobs.py --MEM 9000 --ncores 1 --jobname bamsort_name --command "python ~/software/utils/bam/bamSortName.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .no_duplicates.bam"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_2.txt | python ~/software/utils/submitJobs.py --MEM 9000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortName.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .no_duplicates.bam"
echo "febc_D_ATAC" | python ~/software/utils/submitJobs.py --MEM 18000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortName.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .no_duplicates.bam"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 9000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortName.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .no_duplicates.bam"

#Convert BAMS to bed file of fragments
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bamToFragmentBed --command "python ~/software/utils/bam/bamToFragmentBed.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .sortedByName.bam --outsuffix .fragments.bed.gz --execute True"
cut -f1 failed_fragments.txt | python ~/software/utils/submitJobs.py --MEM 3000 --ncores 2 --jobname bamToFragmentBed --command "python ~/software/utils/bam/bamToFragmentBed.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .sortedByName.bam --outsuffix .fragments.bed.gz --execute True"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 3000 --ncores 2 --jobname bamToFragmentBed --command "python ~/software/utils/bam/bamToFragmentBed.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .sortedByName.bam --outsuffix .fragments.bed.gz --execute True"

#Calculate fragment length distributions for each sample
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt  | python ~/software/utils/submitJobs.py --MEM 1500 --jobname countFragmentLengths --command "python ~/software/utils/coverage/bedCountFragmentLengths.py --indir processed/SL1344 --outdir processed/SL1344"

#Convert BED to bigwig
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bed2bigwig --command "python ~/software/utils/coverage/bed2bigwig.py --indir processed/SL1344/ --outdir processed/SL1344/ --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --execute True"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_2.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bed2bigwig --command "python ~/software/utils/coverage/bed2bigwig.py --indir processed/SL1344/ --outdir processed/SL1344/ --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --execute True"

###Create cutsite coverage plots ####

#Call peaks
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname atacMacsCallPeak --command "python ~/software/utils/coverage/atacMacsCallPeak.py --indir processed/SL1344/ --outdir processed/SL1344/ --execute True"
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names_3.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname atacMacsCallPeak --command "python ~/software/utils/coverage/atacMacsCallPeak.py --indir processed/SL1344/ --outdir processed/SL1344/ --execute True"

#Count the number of reads mapping to consensus peaks
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir processed/SL1344/ --gtf annotations/ATAC_consensus_peaks.gff3 --strand 0 --countsSuffix .consensus_peaks.counts.txt --bamSuffix .sortedByName.bam --execute True --donotsort False --O True"


#### PEER ####
#Run PEER on each condition separately using only expressed genes
echo "PEER" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runPEER --command "python ~/software/utils/runPEER.py --input results/ATAC/PEER/input/naive.exprs.txt --outdir results/ATAC/PEER/output/naive_10/ --n_factors 10"
echo "PEER" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runPEER --command "python ~/software/utils/runPEER.py --input results/ATAC/PEER/input/IFNg.exprs.txt --outdir results/ATAC/PEER/output/IFNg_10/ --n_factors 10"
echo "PEER" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runPEER --command "python ~/software/utils/runPEER.py --input results/ATAC/PEER/input/SL1344.exprs.txt --outdir results/ATAC/PEER/output/SL1344_10/ --n_factors 10"
echo "PEER" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runPEER --command "python ~/software/utils/runPEER.py --input results/ATAC/PEER/input/IFNg_SL1344.exprs.txt --outdir results/ATAC/PEER/output/IFNg_SL1344_10/ --n_factors 10"

#Run PEER on TPM values
echo "PEER" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runPEER --command "python ~/software/utils/runPEER.py --input results/ATAC/PEER/input/naive.exprs_tpm.txt --outdir results/ATAC/PEER/output/naive_tpm_10/ --n_factors 10"
echo "PEER" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runPEER --command "python ~/software/utils/runPEER.py --input results/ATAC/PEER/input/IFNg.exprs_tpm.txt --outdir results/ATAC/PEER/output/IFNg_tpm_10/ --n_factors 10"
echo "PEER" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runPEER --command "python ~/software/utils/runPEER.py --input results/ATAC/PEER/input/SL1344.exprs_tpm.txt --outdir results/ATAC/PEER/output/SL1344_tpm_10/ --n_factors 10"
echo "PEER" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runPEER --command "python ~/software/utils/runPEER.py --input results/ATAC/PEER/input/IFNg_SL1344.exprs_tpm.txt --outdir results/ATAC/PEER/output/IFNg_SL1344_tpm_10/ --n_factors 10"

#### Count allele-specifc expression ####
#Count allele-specific expression
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 10000 --jobname bamCountASE --command "python ~/software/utils/rasqual/bamCountASE.py --indir processed/SL1344/ --outdir processed/SL1344/ --insuffix .no_duplicates.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites ../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True --Xmx 8g"

#Construct a sample-sample map for meregeASECounts.py script
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | awk -v OFS='\t' '{print $1, $1}' > processed/SL1344/sample_sample_map.txt

#Merge all allele-specific counts into one matrix
echo "mergeASECounts" | python ~/software/utils/submitJobs.py --MEM 38000 --queue hugemem --jobname mergeASECounts --command "python ~/software/utils/rasqual/mergeASECounts.py --sample_list processed/SL1344/sample_sample_map.txt --indir processed/SL1344 --suffix .ASEcounts > results/ATAC/ATAC_combined_ASE_counts.txt"

#### Add ASE counts into the vcf files ####
# Extract genotype lists for each condition
cut -f2 results/ATAC/rasqual/input/naive.sg_map.txt > results/ATAC/rasqual/input/naive.genotypes.txt
cut -f2 results/ATAC/rasqual/input/IFNg.sg_map.txt > results/ATAC/rasqual/input/IFNg.genotypes.txt
cut -f2 results/ATAC/rasqual/input/SL1344.sg_map.txt > results/ATAC/rasqual/input/SL1344.genotypes.txt
cut -f2 results/ATAC/rasqual/input/IFNg_SL1344.sg_map.txt > results/ATAC/rasqual/input/IFNg_SL1344.genotypes.txt

#Extract samples from the global VCF file
bcftools view -Oz -S results/ATAC/rasqual/input/naive.genotypes.txt ../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/ATAC/rasqual/input/naive.vcf &
bcftools view -Oz -S results/ATAC/rasqual/input/IFNg.genotypes.txt ../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/ATAC/rasqual/input/IFNg.vcf &
bcftools view -Oz -S results/ATAC/rasqual/input/SL1344.genotypes.txt ../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/ATAC/rasqual/input/SL1344.vcf &
bcftools view -Oz -S results/ATAC/rasqual/input/IFNg_SL1344.genotypes.txt ../macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/ATAC/rasqual/input/IFNg_SL1344.vcf &

#Add ASE counts into the VCF file
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 64000 --jobname vcfAddASE --queue hugemem --command "python ~/software/utils/rasqual/vcfAddASE.py --ASEcounts results/ATAC/ATAC_combined_ASE_counts.txt --ASESampleGenotypeMap results/ATAC/rasqual/input/naive.sg_map.txt --VCFfile results/ATAC/rasqual/input/naive.vcf | bgzip > results/ATAC/rasqual/input/naive.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 64000 --jobname vcfAddASE --queue hugemem --command "python ~/software/utils/rasqual/vcfAddASE.py --ASEcounts results/ATAC/ATAC_combined_ASE_counts.txt --ASESampleGenotypeMap results/ATAC/rasqual/input/IFNg.sg_map.txt --VCFfile results/ATAC/rasqual/input/IFNg.vcf | bgzip > results/ATAC/rasqual/input/IFNg.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 64000 --jobname vcfAddASE --queue hugemem --command "python ~/software/utils/rasqual/vcfAddASE.py --ASEcounts results/ATAC/ATAC_combined_ASE_counts.txt --ASESampleGenotypeMap results/ATAC/rasqual/input/SL1344.sg_map.txt --VCFfile results/ATAC/rasqual/input/SL1344.vcf | bgzip > results/ATAC/rasqual/input/SL1344.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 64000 --jobname vcfAddASE --queue hugemem --command "python ~/software/utils/rasqual/vcfAddASE.py --ASEcounts results/ATAC/ATAC_combined_ASE_counts.txt --ASESampleGenotypeMap results/ATAC/rasqual/input/IFNg_SL1344.sg_map.txt --VCFfile results/ATAC/rasqual/input/IFNg_SL1344.vcf | bgzip > results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz"

#Index the VCF files
tabix -p vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz &
tabix -p vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz &
tabix -p vcf results/ATAC/rasqual/input/SL1344.ASE.vcf.gz &
tabix -p vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz &

#RUN rasqual on all peaks
#naive
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.gc_library_size.bin --n 42 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/naive_100kb/batches/naive_100kb --covariates results/ATAC/rasqual/input/naive.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5'"
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/IFNg.expression.bin --offsets results/ATAC/rasqual/input/IFNg.gc_library_size.bin --n 41 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/IFNg_100kb/batches/IFNg_100kb --covariates results/ATAC/rasqual/input/IFNg.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5'"
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/SL1344.expression.bin --offsets results/ATAC/rasqual/input/SL1344.gc_library_size.bin --n 31 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/SL1344.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/SL1344_100kb/batches/SL1344_100kb --covariates results/ATAC/rasqual/input/SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5'"
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/IFNg_SL1344.expression.bin --offsets results/ATAC/rasqual/input/IFNg_SL1344.gc_library_size.bin --n 31 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches/IFNg_SL1344_100kb --covariates results/ATAC/rasqual/input/IFNg_SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5'"



#Run rasqual using a random permutation
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.gc_library_size.bin --n 42 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/naive_100kb/batches/naive_100kb --covariates results/ATAC/rasqual/input/naive.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5 --random-permutation'"
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/IFNg.expression.bin --offsets results/ATAC/rasqual/input/IFNg.gc_library_size.bin --n 41 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/IFNg_100kb/batches/IFNg_100kb --covariates results/ATAC/rasqual/input/IFNg.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5 --random-permutation'"
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/SL1344.expression.bin --offsets results/ATAC/rasqual/input/SL1344.gc_library_size.bin --n 31 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/SL1344.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/SL1344_100kb/batches/SL1344_100kb --covariates results/ATAC/rasqual/input/SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5 --random-permutation'"
cat results/ATAC/rasqual/input/all_peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/IFNg_SL1344.expression.bin --offsets results/ATAC/rasqual/input/IFNg_SL1344.gc_library_size.bin --n 31 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches/IFNg_SL1344_100kb --covariates results/ATAC/rasqual/input/IFNg_SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5 --random-permutation'"




#Rerun failed peaks
cat results/ATAC/rasqual/input/IFNg_SL1344_failed_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_naive --ncores 5 --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/ATAC/rasqual/input/IFNg_SL1344.expression.bin --offsets results/ATAC/rasqual/input/IFNg_SL1344.gc_library_size.bin --n 31 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_100kb.txt --outprefix results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches/IFNg_SL1344_100kb --covariates results/ATAC/rasqual/input/IFNg_SL1344.covariates.bin --execute True --rasqualBin rasqual --parameters '\--force --n-threads 5'"

#Merge batches into single files
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/ATAC/rasqual/output/naive_100kb/batches/naive_100kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/ATAC/rasqual/output/IFNg_100kb/batches/IFNg_100kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/ATAC/rasqual/output/SL1344_100kb/batches/SL1344_100kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches/IFNg_SL1344_100kb"

#Move merged files to parent directory
mv results/ATAC/rasqual/output/naive_100kb/batches/naive_100kb.txt results/ATAC/rasqual/output/naive_100kb/
mv results/ATAC/rasqual/output/IFNg_100kb/batches/IFNg_100kb.txt results/ATAC/rasqual/output/IFNg_100kb/
mv results/ATAC/rasqual/output/SL1344_100kb/batches/SL1344_100kb.txt results/ATAC/rasqual/output/SL1344_100kb/
mv results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches/IFNg_SL1344_100kb.txt results/ATAC/rasqual/output/IFNg_SL1344_100kb/

# Compress batches
tar czf results/ATAC/rasqual/output/naive_100kb/batches.tar.gz results/ATAC/rasqual/output/naive_100kb/batches/ &
tar czf results/ATAC/rasqual/output/IFNg_100kb/batches.tar.gz results/ATAC/rasqual/output/IFNg_100kb/batches/ &
tar czf results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches.tar.gz results/ATAC/rasqual/output/IFNg_SL1344_100kb/batches/ &
tar czf results/ATAC/rasqual/output/SL1344_100kb/batches.tar.gz results/ATAC/rasqual/output/SL1344_100kb/batches/ &

#Remove uncompressed batches
rm -r results/ATAC/rasqual/output/*_100kb/batches/

#Extract completed gene ids
cut -f1 results/ATAC/rasqual/output/naive_100kb/naive_100kb.txt | uniq > results/ATAC/rasqual/output/naive_100kb/naive_100kb.completed_ids.txt &
cut -f1 results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.txt | uniq > results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.completed_ids.txt &
cut -f1 results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.txt | uniq > results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.completed_ids.txt &
cut -f1 results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.txt | uniq > results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.completed_ids.txt &

#Sort and index rasqual output for tabix
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED  results/ATAC/rasqual/output/naive_100kb/naive_100kb.txt | sort -k3,3 -k4,4n | bgzip > results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.txt | sort -k3,3 -k4,4n | bgzip > results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout  "grep -v SKIPPED  results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.txt | sort -k3,3 -k4,4n | bgzip > results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout  "grep -v SKIPPED  results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.txt | sort -k3,3 -k4,4n | bgzip > results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz"

#Index text files with tabix
tabix -s3 -b4 -e4 -f results/ATAC/rasqual/output/naive_100kb/naive_100kb.sorted.txt.gz &
tabix -s3 -b4 -e4 -f results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.sorted.txt.gz &
tabix -s3 -b4 -e4 -f results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.sorted.txt.gz &
tabix -s3 -b4 -e4 -f results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.sorted.txt.gz &

### eigenMT ####
#Convert rasqual output into format suitable for eigenMT
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/ATAC/rasqual/output/naive_100kb/naive_100kb.txt > results/ATAC/rasqual/output/naive_100kb/naive_100kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.txt > results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.txt > results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.txt > results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.eigenMT_input.txt"

#Remove the uncompressed output files
rm results/ATAC/rasqual/output/*_100kb/*_100kb.txt

#Run eigenMT chromosome-by-chromosme
#Make sure that peak locations are properly exported
cat ../macrophage-gxe-study/macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 4 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir ../macrophage-gxe-study/results/SL1344/eigenMT/input/ --genepos results/ATAC/eigenMT/input/gene_positions.txt --QTL results/ATAC/rasqual/output/naive_100kb/naive_100kb.eigenMT_input.txt --out_prefix results/ATAC/rasqual/output/naive_100kb/naive_50kb --cis_dist 5e4 --eigenMT_path ~/software/eigenMT/eigenMT.py --external"
cat ../macrophage-gxe-study/macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 4 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir ../macrophage-gxe-study/results/SL1344/eigenMT/input/ --genepos results/ATAC/eigenMT/input/gene_positions.txt --QTL results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.eigenMT_input.txt --out_prefix results/ATAC/rasqual/output/IFNg_100kb/IFNg_50kb --cis_dist 5e4 --eigenMT_path ~/software/eigenMT/eigenMT.py --external"
cat ../macrophage-gxe-study/macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 4 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir ../macrophage-gxe-study/results/SL1344/eigenMT/input/ --genepos results/ATAC/eigenMT/input/gene_positions.txt --QTL results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.eigenMT_input.txt --out_prefix results/ATAC/rasqual/output/SL1344_100kb/SL1344_50kb --cis_dist 5e4 --eigenMT_path ~/software/eigenMT/eigenMT.py --external"
cat ../macrophage-gxe-study/macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 4 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --chromosome_dir ../macrophage-gxe-study/results/SL1344/eigenMT/input/ --genepos results/ATAC/eigenMT/input/gene_positions.txt --QTL results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.eigenMT_input.txt --out_prefix results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_50kb --cis_dist 5e4 --eigenMT_path ~/software/eigenMT/eigenMT.py --external"

#Merge eigenMT output 
cat results/ATAC/rasqual/output/naive_100kb/naive_50kb.chr_*.eigenMT.txt | grep -v snps > results/ATAC/rasqual/output/naive_100kb/naive_50kb.eigenMT.txt
cat results/ATAC/rasqual/output/IFNg_100kb/IFNg_50kb.chr_*.eigenMT.txt | grep -v snps > results/ATAC/rasqual/output/IFNg_100kb/IFNg_50kb.eigenMT.txt
cat results/ATAC/rasqual/output/SL1344_100kb/SL1344_50kb.chr_*.eigenMT.txt | grep -v snps > results/ATAC/rasqual/output/SL1344_100kb/SL1344_50kb.eigenMT.txt
cat results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_50kb.chr_*.eigenMT.txt | grep -v snps > results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_50kb.eigenMT.txt

#Compress eigenMT input 
gzip results/ATAC/rasqual/output/naive_100kb/naive_100kb.eigenMT_input.txt &
gzip results/ATAC/rasqual/output/IFNg_100kb/IFNg_100kb.eigenMT_input.txt &
gzip results/ATAC/rasqual/output/SL1344_100kb/SL1344_100kb.eigenMT_input.txt &
gzip results/ATAC/rasqual/output/IFNg_SL1344_100kb/IFNg_SL1344_100kb.eigenMT_input.txt &


##### Run FASTQTL #####
#prepare expression data
bgzip results/ATAC/fastqtl/input/naive.expression_cqn.txt && tabix -p bed results/ATAC/fastqtl/input/naive.expression_cqn.txt.gz
bgzip results/ATAC/fastqtl/input/IFNg.expression_cqn.txt && tabix -p bed results/ATAC/fastqtl/input/IFNg.expression_cqn.txt.gz
bgzip results/ATAC/fastqtl/input/SL1344.expression_cqn.txt && tabix -p bed results/ATAC/fastqtl/input/SL1344.expression_cqn.txt.gz
bgzip results/ATAC/fastqtl/input/IFNg_SL1344.expression_cqn.txt && tabix -p bed results/ATAC/fastqtl/input/IFNg_SL1344.expression_cqn.txt.gz

#Run FastQTL with CQN-normalised data and covariates (3 PCs + sex) (50kb window)
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/naive.covariates_cqn.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output/naive_50kb_cqn_perm --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --bed results/ATAC/fastqtl/input/IFNg.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/IFNg.covariates_cqn.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output/IFNg_50kb_cqn_perm --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/SL1344.ASE.vcf.gz --bed results/ATAC/fastqtl/input/SL1344.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/SL1344.covariates_cqn.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output/SL1344_50kb_cqn_perm --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz --bed results/ATAC/fastqtl/input/IFNg_SL1344.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/IFNg_SL1344.covariates_cqn.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output/IFNg_SL1344_50kb_cqn_perm --execute True"

#Merge batches
zcat results/ATAC/fastqtl/output/naive_50kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/naive_50kb_cqn_perm.txt.gz
zcat results/ATAC/fastqtl/output/IFNg_50kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_50kb_cqn_perm.txt.gz
zcat results/ATAC/fastqtl/output/SL1344_50kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/SL1344_50kb_cqn_perm.txt.gz
zcat results/ATAC/fastqtl/output/IFNg_SL1344_50kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_SL1344_50kb_cqn_perm.txt.gz
rm results/ATAC/fastqtl/output/*chunk_*

#Run FastQTL with CQN-normalised data and covariates (3 PCs + sex) (100kb window)
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/naive.covariates_cqn.txt --W 100000 --permute '100 10000' --out results/ATAC/fastqtl/output/naive_100kb_cqn_perm --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --bed results/ATAC/fastqtl/input/IFNg.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/IFNg.covariates_cqn.txt --W 100000 --permute '100 10000' --out results/ATAC/fastqtl/output/IFNg_100kb_cqn_perm --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/SL1344.ASE.vcf.gz --bed results/ATAC/fastqtl/input/SL1344.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/SL1344.covariates_cqn.txt --W 100000 --permute '100 10000' --out results/ATAC/fastqtl/output/SL1344_100kb_cqn_perm --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz --bed results/ATAC/fastqtl/input/IFNg_SL1344.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/IFNg_SL1344.covariates_cqn.txt --W 100000 --permute '100 10000' --out results/ATAC/fastqtl/output/IFNg_SL1344_100kb_cqn_perm --execute True"

#Merge batches
zcat results/ATAC/fastqtl/output/naive_100kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/naive_100kb_cqn_perm.txt.gz
zcat results/ATAC/fastqtl/output/IFNg_100kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_100kb_cqn_perm.txt.gz
zcat results/ATAC/fastqtl/output/SL1344_100kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/SL1344_100kb_cqn_perm.txt.gz
zcat results/ATAC/fastqtl/output/IFNg_SL1344_100kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_SL1344_100kb_cqn_perm.txt.gz
rm results/ATAC/fastqtl/output/*chunk_*

#Get fastqtl results for all SNP-peak pairs in a 100kb window
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/naive.covariates_cqn.txt --W 100000 --out results/ATAC/fastqtl/output/naive_100kb_full --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --bed results/ATAC/fastqtl/input/IFNg.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/IFNg.covariates_cqn.txt --W 100000 --out results/ATAC/fastqtl/output/IFNg_100kb_full --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/SL1344.ASE.vcf.gz --bed results/ATAC/fastqtl/input/SL1344.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/SL1344.covariates_cqn.txt --W 100000 --out results/ATAC/fastqtl/output/SL1344_100kb_full --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz --bed results/ATAC/fastqtl/input/IFNg_SL1344.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/IFNg_SL1344.covariates_cqn.txt --W 100000 --out results/ATAC/fastqtl/output/IFNg_SL1344_100kb_full --execute True"

#Merge chunks into single files
zcat results/ATAC/fastqtl/output/naive_100kb_full.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/naive_100kb_pvalues.txt.gz &
zcat results/ATAC/fastqtl/output/IFNg_100kb_full.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_100kb_pvalues.txt.gz &
zcat results/ATAC/fastqtl/output/SL1344_100kb_full.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/SL1344_100kb_pvalues.txt.gz &
zcat results/ATAC/fastqtl/output/IFNg_SL1344_100kb_full.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_SL1344_100kb_pvalues.txt.gz &

#Remove chunks
rm results/ATAC/fastqtl/output/*.chunk_*

#Add SNP coordinates
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --fastqtl results/ATAC/fastqtl/output/naive_100kb_pvalues.txt.gz | bgzip > results/ATAC/fastqtl/output/naive_100kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --fastqtl results/ATAC/fastqtl/output/IFNg_100kb_pvalues.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_100kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/ATAC/rasqual/input/SL1344.ASE.vcf.gz --fastqtl results/ATAC/fastqtl/output/SL1344_100kb_pvalues.txt.gz | bgzip > results/ATAC/fastqtl/output/SL1344_100kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz --fastqtl results/ATAC/fastqtl/output/IFNg_SL1344_100kb_pvalues.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_SL1344_100kb_pvalues.coords.txt.gz"

#Sort files by SNP coordinates
#awk command is necessary to change field separator from space to tab
zcat results/ATAC/fastqtl/output/naive_100kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/ATAC/fastqtl/output/naive_100kb_pvalues.sorted.txt.gz
zcat results/ATAC/fastqtl/output/IFNg_100kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/ATAC/fastqtl/output/IFNg_100kb_pvalues.sorted.txt.gz
zcat results/ATAC/fastqtl/output/SL1344_100kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/ATAC/fastqtl/output/SL1344_100kb_pvalues.sorted.txt.gz
zcat results/ATAC/fastqtl/output/IFNg_SL1344_100kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/ATAC/fastqtl/output/IFNg_SL1344_100kb_pvalues.sorted.txt.gz
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 -o FarmOut/sort_fastqtl.%J.jobout "./sort_fastqtl.sh"

#Index the output files using Tabix
tabix -s2 -b3 -e3 -f results/ATAC/fastqtl/output/naive_100kb_pvalues.sorted.txt.gz
tabix -s2 -b3 -e3 -f results/ATAC/fastqtl/output/IFNg_100kb_pvalues.sorted.txt.gz
tabix -s2 -b3 -e3 -f results/ATAC/fastqtl/output/SL1344_100kb_pvalues.sorted.txt.gz
tabix -s2 -b3 -e3 -f results/ATAC/fastqtl/output/IFNg_SL1344_100kb_pvalues.sorted.txt.gz


#Run FastQTL with CQN-normalised data and covariates (3 PCs + sex) (500kb window)

#Get fastqtl results for all SNP-peak pairs in a 500kb window
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/naive.covariates_cqn.txt --W 500000 --out results/ATAC/fastqtl/output/naive_500kb_full --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --bed results/ATAC/fastqtl/input/IFNg.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/IFNg.covariates_cqn.txt --W 500000 --out results/ATAC/fastqtl/output/IFNg_500kb_full --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/SL1344.ASE.vcf.gz --bed results/ATAC/fastqtl/input/SL1344.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/SL1344.covariates_cqn.txt --W 500000 --out results/ATAC/fastqtl/output/SL1344_500kb_full --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz --bed results/ATAC/fastqtl/input/IFNg_SL1344.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/IFNg_SL1344.covariates_cqn.txt --W 500000 --out results/ATAC/fastqtl/output/IFNg_SL1344_500kb_full --execute True"

#Merge chunks into single files
zcat results/ATAC/fastqtl/output/naive_500kb_full.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/naive_500kb_pvalues.txt.gz &
zcat results/ATAC/fastqtl/output/IFNg_500kb_full.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_500kb_pvalues.txt.gz &
zcat results/ATAC/fastqtl/output/SL1344_500kb_full.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/SL1344_500kb_pvalues.txt.gz &
zcat results/ATAC/fastqtl/output/IFNg_SL1344_500kb_full.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_SL1344_500kb_pvalues.txt.gz &

#Remove chunks
rm results/ATAC/fastqtl/output/*.chunk_*

#Add SNP coordinates
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --fastqtl results/ATAC/fastqtl/output/naive_500kb_pvalues.txt.gz | bgzip > results/ATAC/fastqtl/output/naive_500kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --fastqtl results/ATAC/fastqtl/output/IFNg_500kb_pvalues.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_500kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/ATAC/rasqual/input/SL1344.ASE.vcf.gz --fastqtl results/ATAC/fastqtl/output/SL1344_500kb_pvalues.txt.gz | bgzip > results/ATAC/fastqtl/output/SL1344_500kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.vcf.gz --fastqtl results/ATAC/fastqtl/output/IFNg_SL1344_500kb_pvalues.txt.gz | bgzip > results/ATAC/fastqtl/output/IFNg_SL1344_500kb_pvalues.coords.txt.gz"

#Sort files by SNP coordinates
#awk command is necessary to change field separator from space to tab
zcat results/ATAC/fastqtl/output/naive_500kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/ATAC/fastqtl/output/naive_500kb_pvalues.sorted.txt.gz
zcat results/ATAC/fastqtl/output/IFNg_500kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/ATAC/fastqtl/output/IFNg_500kb_pvalues.sorted.txt.gz
zcat results/ATAC/fastqtl/output/SL1344_500kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/ATAC/fastqtl/output/SL1344_500kb_pvalues.sorted.txt.gz
zcat results/ATAC/fastqtl/output/IFNg_SL1344_500kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/ATAC/fastqtl/output/IFNg_SL1344_500kb_pvalues.sorted.txt.gz
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 -o FarmOut/sort_fastqtl.%J.jobout "./sort_fastqtl.sh"

#Index the output files using Tabix
tabix -s2 -b3 -e3 -f results/ATAC/fastqtl/output/naive_500kb_pvalues.sorted.txt.gz
tabix -s2 -b3 -e3 -f results/ATAC/fastqtl/output/IFNg_500kb_pvalues.sorted.txt.gz
tabix -s2 -b3 -e3 -f results/ATAC/fastqtl/output/SL1344_500kb_pvalues.sorted.txt.gz
tabix -s2 -b3 -e3 -f results/ATAC/fastqtl/output/IFNg_SL1344_500kb_pvalues.sorted.txt.gz


#Extract FASTA sequence corresponding to peaks
/software/CGP/external-apps/cufflinks-1.3.0/bin/gffread ATAC_consensus_peaks.gff3  -g ../../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -w ATAC_consnesus_peaks.fasta

#Calculate nucleotide content for each peak
bedtools nuc -fi ../../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ATAC_consensus_peaks.gff3  > ATAC_consensus_peaks.nuc_content.txt

##### FASTQTL single permutation run #####
bcftools query -l naive.ASE.vcf.gz | sort -R > naive.sample_names.txt
bcftools reheader -s naive.sample_names.txt naive.ASE.vcf.gz > naive.ASE.permuted.vcf.gz
tabix -p vcf naive.ASE.permuted.vcf.gz

bcftools query -l IFNg.ASE.vcf.gz | sort -R > IFNg.sample_names.txt
bcftools reheader -s IFNg.sample_names.txt IFNg.ASE.vcf.gz > IFNg.ASE.permuted.vcf.gz
tabix -p vcf IFNg.ASE.permuted.vcf.gz

bcftools query -l SL1344.ASE.vcf.gz | sort -R > SL1344.sample_names.txt
bcftools reheader -s SL1344.sample_names.txt SL1344.ASE.vcf.gz > SL1344.ASE.permuted.vcf.gz
tabix -p vcf SL1344.ASE.permuted.vcf.gz

bcftools query -l IFNg_SL1344.ASE.vcf.gz | sort -R > IFNg_SL1344.sample_names.txt
bcftools reheader -s IFNg_SL1344.sample_names.txt IFNg_SL1344.ASE.vcf.gz > IFNg_SL1344.ASE.permuted.vcf.gz
tabix -p vcf IFNg_SL1344.ASE.permuted.vcf.gz

#Run FastQTL with CQN-normalised data and covariates (3 PCs + sex) (50kb window)
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.permuted.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/naive.covariates_cqn.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output_permutation/naive_50kb_cqn_perm --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/IFNg.ASE.permuted.vcf.gz --bed results/ATAC/fastqtl/input/IFNg.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/IFNg.covariates_cqn.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output_permutation/IFNg_50kb_cqn_perm --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/SL1344.ASE.permuted.vcf.gz --bed results/ATAC/fastqtl/input/SL1344.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/SL1344.covariates_cqn.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output_permutation/SL1344_50kb_cqn_perm --execute True"
cat results/ATAC/fastqtl/input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/IFNg_SL1344.ASE.permuted.vcf.gz --bed results/ATAC/fastqtl/input/IFNg_SL1344.expression_cqn.txt.gz --cov results/ATAC/fastqtl/input/IFNg_SL1344.covariates_cqn.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output_permutation/IFNg_SL1344_50kb_cqn_perm --execute True"

#Merge batches
zcat results/ATAC/fastqtl/output_permutation/naive_50kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output_permutation/naive_50kb_cqn_perm.txt.gz
zcat results/ATAC/fastqtl/output_permutation/IFNg_50kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output_permutation/IFNg_50kb_cqn_perm.txt.gz
zcat results/ATAC/fastqtl/output_permutation/SL1344_50kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output_permutation/SL1344_50kb_cqn_perm.txt.gz
zcat results/ATAC/fastqtl/output_permutation/IFNg_SL1344_50kb_cqn_perm.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output_permutation/IFNg_SL1344_50kb_cqn_perm.txt.gz
rm results/ATAC/fastqtl/output_permutation/*chunk_*


##### MOTIFS #####
#Scan all ATAC peaks using FIMO
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q long -M 1000 -o FarmOut/FIMO_scan_long.%J.jobout "~/software/meme/bin/fimo --max-stored-scores 100000000 --thresh 1e-5 --text --verbosity 2 annotations/motif_databases/CIS-BP/Homo_sapiens.meme annotations/ATAC_consensus_peaks.fasta > results/ATAC/FIMO_CISBP_results.long.txt" 
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q long -M 1000 -o FarmOut/FIMO_scan_long.%J.jobout "~/software/meme/bin/fimo --max-stored-scores 100000000 --thresh 1e-4 --text --verbosity 2 annotations/motif_databases/CIS-BP/Homo_sapiens.meme annotations/ATAC_consensus_peaks.fasta > results/ATAC/FIMO_CISBP_results.1e-4.txt" 

#Scan promoter regions as a background
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q long -M 1000 -o FarmOut/FIMO_scan_long.%J.jobout "~/software/meme/bin/fimo --max-stored-scores 100000000 --thresh 1e-5 --text --verbosity 2 annotations/motif_databases/CIS-BP/Homo_sapiens.meme annotations/PWMEnrich_promoter_regions.fasta > results/ATAC/FIMO_CISBP_promoter_matches.txt" 

#Scan QTL peaks for motif disruption events
cat motif_batches.txt | python ~/software/utils/submitJobs.py --MEM 4000 --jobname calculate_motif_disruptions --command "/software/R-3.2.2/bin/Rscript macrophage-chromatin/ATAC/finemapping/caQTL_motif_disruptions.R"
cat results/ATAC/motif_analysis/motif_disruption_batch_*.txt > results/ATAC/motif_analysis/motif_disruption.txt
rm results/ATAC/motif_analysis/motif_disruption_batch_*.txt

##### ChIP-seq overlaps #####
#H3K27Ac
gat-run.py --segments=results/ATAC/DA/ATAC_clustered_peaks.bed --annotations=results/public_chromatin/joint_peaks/H3K27Ac_clustered_peaks.bed  --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=1000 --log=results/public_chromatin/annotation_overlaps/H3K27Ac_overlap.gat_log.txt --with-segment-tracks > results/public_chromatin/annotation_overlaps/H3K27Ac_overlap.gat.txt
#STAT1
gat-run.py --segments=results/ATAC/DA/ATAC_clustered_peaks.bed --annotations=results/Ivashkiv/DA/STAT1_grouped_peaks.bed  --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=1000 --log=results/annotation_overlaps/STAT1_overlap.gat_log.txt --with-segment-tracks > results/annotation_overlaps/STAT1_overlap.gat.txt
#IRF1
gat-run.py --segments=results/ATAC/DA/ATAC_clustered_peaks.bed --annotations=results/Ivashkiv/peak_calls/IRF1_joint_peaks.bed --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=1000 --log=results/annotation_overlaps/IRF1_overlap.gat_log.txt --with-segment-tracks > results/annotation_overlaps/IRF1_overlap.gat.txt
#CIITA-RFX5
gat-run.py --segments=results/ATAC/DA/ATAC_clustered_peaks.bed --annotations=results/public_chromatin/joint_peaks/CIITA-RFX5_joint_peaks.bed --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=1000 --log=results/public_chromatin/annotation_overlaps/CIITA-RFX5_overlap.gat_log.txt --with-segment-tracks > results/public_chromatin/annotation_overlaps/CIITA-RFX5_overlap.gat.txt

#Enrichments for PU.1, CEBPb and CTCF
gat-run.py --segments=annotations/ATAC_consensus_peaks.bed  --annotations=results/ATAC/ChIP_enrichment/naive_combined_peaks.bed --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=100 --log=results/ATAC/ChIP_enrichment/gat_output/ATAC_naive_overlap.gat_log.txt --with-segment-tracks > results/ATAC/ChIP_enrichment/gat_output/ATAC_naive_overlap.gat.txt
gat-run.py --segments=results/ATAC/ChIP_enrichment/naive_combined_peaks.bed  --annotations=results/ATAC/ChIP_enrichment/naive_combined_peaks.bed --workspace=../../annotations/blacklists/GRCh38_filtered_gapped_genome.bed --num-samples=100 --log=results/ATAC/ChIP_enrichment/gat_output/naive_naive_overlap.gat_log.txt --with-segment-tracks > results/ATAC/ChIP_enrichment/gat_output/naive_naive_overlap.gat.txt

#### Construct credible sets #####
echo "hello" | python ~/software/utils/submitJobs.py --MEM 14000 --jobname construct_credible_sets --queue long --command "/software/R-3.3.0/bin/Rscript macrophage-gxe-study/ATAC/finemapping/ATAC_construct_credible_sets.R"

