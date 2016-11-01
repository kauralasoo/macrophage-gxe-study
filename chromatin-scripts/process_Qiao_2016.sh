#Convert Ivashkiv data to fastq
cut -f2 macrophage-gxe-study/data/chromatin/ChIP/Qiao_2016_sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname sra2fastq --command "python ~/software/utils/fastq/sra2fastq.py --indir data/Qiao_2016/sra/ --outdir data/Qiao_2016/fastq/"

#Align Chip-seq data to the genome using bwa aln
cut -f1 macrophage-gxe-study/data/chromatin/ChIP/Qiao_2016_sample_names.txt | head -n 4 | python ~/software/utils/submitJobs.py --MEM 10000 --ncores 8 --jobname bwa_aln --command  "python ~/software/utils/align/bwaAlnSE.py --outputDir processed/Qiao_2016/ --fastqDir data/Qiao_2016/fastq/ --nCores 8 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Convert to bam
cut -f1 macrophage-gxe-study/data/chromatin/ChIP/Qiao_2016_sample_names.txt | head -n 4 | python ~/software/utils/submitJobs.py --MEM 10000 --jobname bwa_samse --command "python ~/software/utils/align/bwaSamse.py --outputDir processed/Qiao_2016/ --fastqDir data/Qiao_2016/fastq/ --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Sort bams by coordinate
cut -f1 macrophage-gxe-study/data/chromatin/ChIP/Qiao_2016_sample_names.txt | head -n 4 | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/Qiao_2016/ --outdir processed/Qiao_2016/ --insuffix .aligned.bam"

#Remove BWA header because it clashes with the Picard MarkDuplicates
cut -f1 macrophage-gxe-study/data/chromatin/ChIP/Qiao_2016_sample_names.txt | head -n 4 | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 1 --jobname remove_header --command "python ~/software/utils/bam/bamRemoveBwaHeader.py --indir processed/Qiao_2016/ --outdir processed/Qiao_2016/ --insuffix .sortedByCoord.bam --outsuffix .reheadered.bam"

#Remove duplicates
cut -f1 macrophage-gxe-study/data/chromatin/ChIP/Qiao_2016_sample_names.txt | head -n 4 | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 1 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir  processed/Qiao_2016/ --outdir  processed/Qiao_2016/ --execute True --insuffix .reheadered.bam"

#Index bams
cut -f1 macrophage-gxe-study/data/chromatin/ChIP/Qiao_2016_sample_names.txt | head -n 4 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/Qiao_2016/ --insuffix .no_duplicates.bam --execute True"
