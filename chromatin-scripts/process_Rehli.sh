#Align to the genome using bwa aln
cut -f1 macrophage-chromatin/data/Rehli/Rehli_sample_names.txt |  python ~/software/utils/submitJobs.py --MEM 8000 --ncores 8 --jobname bwa_aln --command  "python ~/software/utils/align/bwaAlnSE.py --outputDir processed/Rehli/ --fastqDir data/Rehli/fastq/ --nCores 8 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Convert to bam
cut -f1 macrophage-chromatin/data/Rehli/Rehli_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 8000 --jobname bwa_samse --command "python ~/software/utils/align/bwaSamse.py --outputDir processed/Rehli/ --fastqDir data/Rehli/fastq/ --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Sort bams by coordinate
cut -f1 macrophage-chromatin/data/Rehli/Rehli_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/Rehli/ --outdir processed/Rehli/ --insuffix .aligned.bam"

#Remove BWA header because it clashes with the Picard MarkDuplicates
cut -f1 macrophage-chromatin/data/Rehli/Rehli_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 1 --jobname remove_header --command "python ~/software/utils/bam/bamRemoveBwaHeader.py --indir processed/Rehli/ --outdir processed/Rehli/ --insuffix .sortedByCoord.bam --outsuffix .reheadered.bam"

#Remove duplicates
cut -f1 macrophage-chromatin/data/Rehli/Rehli_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 1 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir  processed/Rehli/ --outdir  processed/Rehli/ --execute True --insuffix .reheadered.bam"

#Index bams
cut -f1 macrophage-chromatin/data/Rehli/Rehli_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/Rehli/ --insuffix .no_duplicates.bam --execute True"

#Convert BAM to BigWig
cut -f1 macrophage-chromatin/data/Rehli/Rehli_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bamToBigWig --command "python ~/software/utils/coverage/bamToBigWig.py --indir processed/Rehli/ --outdir processed/Rehli/ --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --execute True"

#Clean up
rm processed/Rehli/*/*.aligned.bam
rm processed/Rehli/*/*.sortedByCoord.bam
rm processed/Rehli/*/*.reheadered.bam

#Call narrow and broad peaks
cut -f1 macrophage-chromatin/data/Rehli/Rehli_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakNarrow --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Rehli/ --outdir processed/Rehli/ --execute True --control processed/Rehli/IgG/IgG.no_duplicates.bam --broad False"

cut -f1 macrophage-chromatin/data/Rehli/Rehli_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakBroad --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Rehli/ --outdir processed/Rehli/ --execute True --control processed/Rehli/IgG/IgG.no_duplicates.bam --broad True"
