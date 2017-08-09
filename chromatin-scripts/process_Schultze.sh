#Convert Schultze data to fastq
cut -f1 macrophage-chromatin/data/Schultze/Schultze_samples.txt | tail -n 55 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname sra2fastq --command "python ~/software/utils/fastq/sra2fastq.py --indir data/Schultze_2016/sra/ --outdir data/Schultze_2016/fastq/"
echo 'SRR2016717' | python ~/software/utils/submitJobs.py --MEM 1000 --jobname sra2fastq --command "python ~/software/utils/fastq/sra2fastq.py --indir data/Schultze_2016/sra/ --outdir data/Schultze_2016/fastq/"

#Merge split fastq into joint ones and rename
cat macrophage-chromatin/data/Schultze/Schultze_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_fastq --command "python ~/software/utils/fastq/merge-fastq.py --indir data/Schultze_2016/fastq/ --outdir data/Schultze_2016/fastq/ --suffix .fastq.gz"

#Align to the genome using bwa aln
cut -f1 macrophage-chromatin/data/Schultze/Schultze_sample_names.txt |  python ~/software/utils/submitJobs.py --MEM 18000 --ncores 8 --jobname bwa_aln --command  "python ~/software/utils/align/bwaAlnSE.py --outputDir processed/Schultze_2016/ --fastqDir data/Schultze_2016/fastq/ --nCores 8 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Convert to bam
cut -f1  macrophage-chromatin/data/Schultze/Schultze_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 8000 --jobname bwa_samse --command "python ~/software/utils/align/bwaSamse.py --outputDir processed/Schultze_2016/ --fastqDir data/Schultze_2016/fastq/ --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Sort bams by coordinate
cut -f1  macrophage-chromatin/data/Schultze/Schultze_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/Schultze_2016/ --outdir processed/Schultze_2016/ --insuffix .aligned.bam"

#Remove BWA header because it clashes with the Picard MarkDuplicates
cut -f1 macrophage-chromatin/data/Schultze/Schultze_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 1 --jobname remove_header --command "python ~/software/utils/bam/bamRemoveBwaHeader.py --indir processed/Schultze_2016/ --outdir processed/Schultze_2016/ --insuffix .sortedByCoord.bam --outsuffix .reheadered.bam"

#Remove duplicates
cut -f1 macrophage-chromatin/data/Schultze/Schultze_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 1 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir  processed/Schultze_2016/ --outdir  processed/Schultze_2016/ --execute True --insuffix .reheadered.bam"

#Index bams
cut -f1 macrophage-chromatin/data/Schultze/Schultze_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/Schultze_2016/ --insuffix .no_duplicates.bam --execute True"

#Convert BAM to BigWig
cut -f1 macrophage-chromatin/data/Schultze/Schultze_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bamToBigWig --command "python ~/software/utils/coverage/bamToBigWig.py --indir processed/Schultze_2016/ --outdir processed/Schultze_2016/ --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --execute True"

#Clean up
rm processed/Rehli_2015/*/*.aligned.bam
rm processed/Rehli_2015/*/*.sortedByCoord.bam
rm processed/Rehli_2015/*/*.reheadered.bam

#Call narrow and broad peaks
cut -f1 macrophage-chromatin/data/Schultze/Schultze_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakNarrow --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Schultze_2016/ --outdir processed/Schultze_2016/ --execute True --control processed/Schultze_2016/Input_naive/Input_naive.no_duplicates.bam --broad False"

cut -f1 macrophage-chromatin/data/Schultze/Schultze_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakBroad --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Schultze_2016/ --outdir processed/Schultze_2016/ --execute True --control processed/Schultze_2016/Input_naive/Input_naive.no_duplicates.bam --broad True"

#Count number of reads overlapping ATAC peaks
#Count the number of reads mapping to consensus peaks
cut -f1 macrophage-gxe-study/data/chromatin/ChIP/Schultze_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir processed/Schultze_2016/ --gtf annotations/chromatin/ATAC_consensus_peaks.gff3 --strand 0 --countsSuffix .ATAC_peaks.counts.txt --bamSuffix .no_duplicates.bam --execute True --donotsort False --O True"


