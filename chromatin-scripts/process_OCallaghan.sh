#Merge split fastq into joint ones and rename
cat macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_fastq --command "python ~/software/utils/merge-fastq.py --indir data/OCallaghan_2015/fastq/ --outdir data/OCallaghan_2015/fastq/ --suffix _1.fastq.gz"
cat macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_fastq --command "python ~/software/utils/merge-fastq.py --indir data/OCallaghan_2015/fastq/ --outdir data/OCallaghan_2015/fastq/ --suffix _2.fastq.gz"

### ALIGNMENT AND FILTERING ###
#Use bwa mam to align the reads to the genome (although the dataset has 50bp PE reads)
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 10000 --ncores 8 --jobname bwa_mem --command  "python ~/software/utils/align/bwaAlignPE.py --outputDir processed/OCallaghan/ --fastqDir data/OCallaghan_2015/fastq/ --nCores 8 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Sort bams by coordinate
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/OCallaghan/ --outdir processed/OCallaghan/"

#Index bams
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/OCallaghan/ --insuffix .sortedByCoord.bam --execute True"

#Keep only properly paired reads on primary chromosomes (remove MT and contigs)
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 1 --jobname remove_contigs --command "python ~/software/utils/bam/bamRemoveContigsMT.py --indir processed/OCallaghan/ --outdir processed/OCallaghan/"

#Remove BWA header because it clashes with the Picard MarkDuplicates
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 1 --jobname remove_header --command "python ~/software/utils/bam/bamRemoveBwaHeader.py --indir processed/OCallaghan/ --outdir processed/OCallaghan/"

#Remove duplicates
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 1 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir processed/OCallaghan/ --outdir processed/OCallaghan/ --execute True --insuffix .reheadered.bam"

#Index bams
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/OCallaghan/ --insuffix .no_duplicates.bam --execute True"

#### Create fragment coverage files ####

#Sort bams by name
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 9000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortName.py --indir processed/OCallaghan/ --outdir processed/OCallaghan/ --insuffix .no_duplicates.bam"

#Convert BAMS to bed file of fragments
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bamToFragmentBed --command "python ~/software/utils/bam/bamToFragmentBed.py --indir processed/OCallaghan/ --outdir processed/OCallaghan/ --insuffix .sortedByName.bam --outsuffix .fragments.bed.gz --execute True"

#Convert BED to bigwig
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bed2bigwig --command "python ~/software/utils/coverage/bed2bigwig.py --indir processed/OCallaghan/ --outdir processed/OCallaghan/ --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --execute True"

#Clean-up
rm processed/OCallaghan/*/*.filtered.bam
rm processed/OCallaghan/*/*_2??.bam
rm processed/OCallaghan/*/*.sortedByCoord.bam*
rm processed/OCallaghan/*/*.reheadered.bam
rm processed/OCallaghan/*/*.new_header.txt

#Call narrow and broad peaks
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakNarrow --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/OCallaghan/ --outdir processed/OCallaghan/ --execute True --control processed/OCallaghan/Input_ctrl/Input_ctrl.no_duplicates.bam --broad False"
cut -f1 macrophage-chromatin/data/OCallaghan/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakBroad --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/OCallaghan/ --outdir processed/OCallaghan/ --execute True --control processed/OCallaghan/Input_ctrl/Input_ctrl.no_duplicates.bam --broad True"

#Count reads overlapping ATAC peaks
cut -f1 macrophage-gxe-study/data/chromatin/ChIP/OCallaghan_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir processed/OCallaghan/ --gtf annotations/chromatin/ATAC_consensus_peaks.gff3 --strand 0 --countsSuffix .ATAC_peaks.counts.txt --bamSuffix .no_duplicates.bam --execute True --donotsort False --O True"




