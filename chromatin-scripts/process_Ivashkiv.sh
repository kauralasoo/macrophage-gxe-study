#Build a bwa index of the reference
bsub -G team170 -n1 -R "span[hosts=1] select[mem>8000] rusage[mem=8000]" -M 8000 -o bwa_index.%J.jobout "~/software/bin/bwa index -p bwa_index/GRCh38 dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

#Convert Ivashkiv data to fastq
cut -f2 macrophage-chromatin/data/Ivashkiv_2013.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname sra2fastq --command "python ~/software/utils/sra2fastq.py --indir data/Ivashkiv_2013/sra/ --outdir data/Ivashkiv_2013/fastq/"

#Align to the genome using bwa aln
cat macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | tail -n20 | python ~/software/utils/submitJobs.py --MEM 10000 --ncores 8 --jobname bwa_aln --command  "python ~/software/utils/align/bwaAlnSE.py --outputDir processed/Ivashkiv/ --fastqDir data/Ivashkiv_2013/fastq/ --nCores 8 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"
echo "STAT1_rep1_D" | python ~/software/utils/submitJobs.py --MEM 10000 --ncores 8 --jobname bwa_aln --command  "python ~/software/utils/align/bwaAlnSE.py --outputDir processed/Ivashkiv/ --fastqDir data/Ivashkiv_2013/fastq/ --nCores 8 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Convert to bam
cat macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | tail -n20 | python ~/software/utils/submitJobs.py --MEM 10000 --jobname bwa_samse --command "python ~/software/utils/align/bwaSamse.py --outputDir processed/Ivashkiv/ --fastqDir data/Ivashkiv_2013/fastq/ --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Sort bams by coordinate
cat macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | tail -n20 | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/Ivashkiv/ --outdir processed/Ivashkiv/ --insuffix .aligned.bam"

#Remove BWA header because it clashes with the Picard MarkDuplicates
cat macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | tail -n20 | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 1 --jobname remove_header --command "python ~/software/utils/bam/bamRemoveBwaHeader.py --indir processed/Ivashkiv/ --outdir processed/Ivashkiv/ --insuffix .sortedByCoord.bam --outsuffix .reheadered.bam"

#Remove duplicates
cat macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | tail -n20 | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 1 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir  processed/Ivashkiv/ --outdir  processed/Ivashkiv/ --execute True --insuffix .reheadered.bam"

#Index bams
cat macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/Ivashkiv/ --insuffix .no_duplicates.bam --execute True"

#Convert BAM to BigWig
cat macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | head -n1 | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bamToBigWig --command "python ~/software/utils/coverage/bamToBigWig.py --indir processed/Ivashkiv/ --outdir processed/Ivashkiv/ --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --execute True"

#Clean up
rm processed/Ivashkiv/*/*.aligned.bam
rm processed/Ivashkiv/*/*.sortedByCoord.bam
rm processed/Ivashkiv/*/*.reheadered.bam

#Call broad peaks
cut -f1 macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | grep H3K27Ac | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeak --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Ivashkiv/ --outdir processed/Ivashkiv/ --execute True --control processed/Ivashkiv/H3K27ac_Input/H3K27ac_Input.no_duplicates.bam --broad True"

cut -f1 macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | grep STAT1_rep | python ~/software/utils/submitJobs.py --MEM 5000 --jobname chipMacsCallPeak --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Ivashkiv/ --outdir processed/Ivashkiv/ --execute True --control processed/Ivashkiv/STAT1_Input/STAT1_Input.no_duplicates.bam --broad True"

cut -f1 macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | grep IRF1 | grep -v IRF1_Input | python ~/software/utils/submitJobs.py --MEM 5000 --jobname chipMacsCallPeak --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Ivashkiv/ --outdir processed/Ivashkiv/ --execute True --control processed/Ivashkiv/IRF1_Input/IRF1_Input.no_duplicates.bam --broad True"

#Call narrow peaks
cut -f1 macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | grep H3K27Ac | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakNarrow --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Ivashkiv/ --outdir processed/Ivashkiv/ --execute True --control processed/Ivashkiv/H3K27ac_Input/H3K27ac_Input.no_duplicates.bam --broad False"

cut -f1 macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | grep STAT1_rep | python ~/software/utils/submitJobs.py --MEM 5000 --jobname chipMacsCallPeak --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Ivashkiv/ --outdir processed/Ivashkiv/ --execute True --control processed/Ivashkiv/STAT1_Input/STAT1_Input.no_duplicates.bam --broad False"

cut -f1 macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | grep IRF1 | grep -v IRF1_Input | python ~/software/utils/submitJobs.py --MEM 5000 --jobname chipMacsCallPeak --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Ivashkiv/ --outdir processed/Ivashkiv/ --execute True --control processed/Ivashkiv/IRF1_Input/IRF1_Input.no_duplicates.bam --broad False"

#Count the number of reads overlapping joint peak calls
cut -f1 macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | grep H3K27Ac  | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir processed/Ivashkiv/ --gtf annotations/H3K27Ac_joint_peaks.gff3 --strand 0 --countsSuffix .H3K27Ac_joint_peaks.counts.txt --bamSuffix .no_duplicates.bam --execute True --unpaired True"

cut -f1 macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | grep IRF1 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir processed/Ivashkiv/ --gtf annotations/IRF1_joint_peaks.gff3 --strand 0 --countsSuffix .IRF1_joint_peaks.counts.txt --bamSuffix .no_duplicates.bam --execute True --unpaired True"

cut -f1 macrophage-chromatin/data/Ivashkiv/Ivashkiv_sample_names.txt | grep STAT1_rep | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir processed/Ivashkiv/ --gtf results/Ivashkiv/peak_calls/STAT1_joint_peaks.gff3  --strand 0 --countsSuffix .STAT1_joint_peaks.counts.txt --bamSuffix .no_duplicates.bam --execute True --unpaired True"



