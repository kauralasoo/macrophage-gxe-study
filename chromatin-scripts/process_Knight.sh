
### ALIGNMENT AND FILTERING ###
#Use bwa mam to align the reads to the genome (although the dataset has 51bp PE reads)
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 18000  --ncores 3 --jobname bwa_mem --command  "python ~/software/utils/align/bwaAlignPE.py --outputDir processed/Knight/ --fastqDir data/Knight_2014/fastq/ --nCores 2 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Realign failed
echo "Bcell_rep2_naive_input" | python ~/software/utils/submitJobs.py --MEM 27000  --ncores 5 --jobname bwa_mem --command  "python ~/software/utils/align/bwaAlignPE.py --outputDir processed/Knight/ --fastqDir data/Knight_2014/fastq/ --nCores 4 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"
echo "MO_rep1_naive_input" | python ~/software/utils/submitJobs.py --MEM 27000  --ncores 5 --jobname bwa_mem --command  "python ~/software/utils/align/bwaAlignPE.py --outputDir processed/Knight/ --fastqDir data/Knight_2014/fastq/ --nCores 4 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"
echo "MO_rep2_naive_input" | python ~/software/utils/submitJobs.py --MEM 27000  --ncores 5 --jobname bwa_mem --command  "python ~/software/utils/align/bwaAlignPE.py --outputDir processed/Knight/ --fastqDir data/Knight_2014/fastq/ --nCores 4 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Sort bams by coordinate
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/Knight/ --outdir processed/Knight/"

#Index bams
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt| python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/Knight/ --insuffix .sortedByCoord.bam --execute True"

#Keep only properly paired reads on primary chromosomes (remove MT and contigs)
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 1 --jobname remove_contigs --command "python ~/software/utils/bam/bamRemoveContigsMT.py --indir processed/Knight/ --outdir processed/Knight/"

#Remove BWA header because it clashes with the Picard MarkDuplicates
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 1 --jobname remove_header --command "python ~/software/utils/bam/bamRemoveBwaHeader.py --indir processed/Knight/ --outdir processed/Knight/"

#Remove duplicates
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 1 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir processed/Knight/ --outdir processed/Knight/ --execute True --insuffix .reheadered.bam"

#Index bams
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/Knight/ --insuffix .no_duplicates.bam --execute True"

#### Create fragment coverage files ####

#Sort bams by name
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 9000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortName.py --indir processed/Knight/ --outdir processed/Knight/ --insuffix .no_duplicates.bam"

#Convert BAMS to bed file of fragments
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bamToFragmentBed --command "python ~/software/utils/bam/bamToFragmentBed.py --indir processed/Knight/ --outdir processed/Knight/ --insuffix .sortedByName.bam --outsuffix .fragments.bed.gz --execute True"

#Convert BED to bigwig
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bed2bigwig --command "python ~/software/utils/coverage/bed2bigwig.py --indir processed/Knight/ --outdir processed/Knight/ --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --execute True"

#Clean-up
rm processed/Knight/*/*.filtered.bam
rm processed/OCallaghan/*/*_2??.bam
rm processed/Knight/*/*.sortedByCoord.bam*
rm processed/Knight/*/*.reheadered.bam
rm processed/Knight/*/*.new_header.txt


#Call narrow and broad peaks
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakNarrow --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Knight/ --outdir processed/Knight/ --execute True --control processed/Knight/MO_rep1_naive_input/MO_rep1_naive_input.no_duplicates.bam --broad False"
cut -f1 macrophage-chromatin/data/Knight/Knight_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakBroad --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Knight/ --outdir processed/Knight/ --execute True --control processed/Knight/MO_rep1_naive_input/MO_rep1_naive_input.no_duplicates.bam --broad True"
