#Trim adapters using Skewer
echo 'bima_A_ATAC' | python ~/software/utils/submitJobs.py --MEM 1000 --jobname fastq_trim --command "python ~/software/utils/fastq/ATAC_trim_adapters.py --fastqDirIn data/ATAC_Salmonella/fastq --fastqDirOut data/ATAC_Salmonella/fastq_trimmed/ --atacPrimers macrophage-chromatin/data/SL1344/ATAC_index_primers.txt --sampleIndexMap macrophage-chromatin/data/SL1344/ATAC_sample_indexes.txt"

#Align reads with 
echo "bima_A_ATAC" | python ~/software/utils/submitJobs.py --MEM 10000 --ncores 8 --jobname bwa_mem --command  "python ~/software/utils/align/bwaAlignPE.py --outputDir processed/trim_test/ --fastqDir data/ATAC_Salmonella/fastq_trimmed/ --nCores 8 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38 --fastqSuffix .trimmed.fastq.gz"
#Sort by coordinates
echo "bima_A_ATAC" | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/trim_test/ --outdir processed/trim_test/"

#Index bams
echo "bima_A_ATAC"  | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/trim_test/ --insuffix .sortedByCoord.bam --execute True"

#Keep only properly paired reads on primary chromosomes (remove MT and contigs)
echo "bima_A_ATAC"  | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 1 --jobname remove_contigs --command "python ~/software/utils/bam/bamRemoveContigsMT.py --indir processed/trim_test/ --outdir processed/trim_test/"

#Remove BWA header because it clashes with the Picard MarkDuplicates
echo "bima_A_ATAC"  | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 1 --jobname remove_header --command "python ~/software/utils/bam/bamRemoveBwaHeader.py --indir processed/trim_test/ --outdir processed/trim_test/"

#Remove MarkDuplicates
echo "bima_A_ATAC" | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 1 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir processed/trim_test/ --outdir processed/trim_test/ --execute True --insuffix .reheadered.bam"

#Sort by name again
echo "bima_A_ATAC" | python ~/software/utils/submitJobs.py --MEM 18000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortName.py --indir processed/trim_test/ --outdir processed/trim_test/ --insuffix .no_duplicates.bam"

#Convert BAM to BED of fragments
echo "bima_A_ATAC"  | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bamToFragmentBed --command "python ~/software/utils/bam/bamToFragmentBed.py --indir processed/trim_test/ --outdir processed/trim_test/ --insuffix .sortedByName.bam --outsuffix .fragments.bed.gz --execute True"

#Calculate fragment length histogram
zcat bima_A_ATAC.fragments.bed.gz | cut -f5| sort | uniq -c > bima_A_ATAC.fragment_lengts.txt
