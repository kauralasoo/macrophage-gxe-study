#Convert BAMs to FASTQ
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname FarifaxBamToFastq --command "python ~/software/utils/bam/cramToFastq.py --inputDir data/Fairfax/bam/ --outputDir data/Fairfax/fastq/ --inputformat bam"

#Use bwa mem to align the reads to the genome (although the dataset has 50bp PE reads)
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 26000 --ncores 4 --jobname bwa_mem --command  "python ~/software/utils/align/bwaAlignPE.py --outputDir processed/Fairfax/ --fastqDir data/Fairfax/fastq/ --nCores 4 --genomeDir ../../annotations/GRCh38/bwa_index/GRCh38"

#Sort bams by coordinate
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 4000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortCoord.py --indir processed/Fairfax/ --outdir processed/Fairfax/"

#Index bams
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/Fairfax/ --insuffix .sortedByCoord.bam --execute True"

#Keep only properly paired reads on primary chromosomes (remove MT and contigs)
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 2000 --ncores 1 --jobname remove_contigs --command "python ~/software/utils/bam/bamRemoveContigsMT.py --indir processed/Fairfax/ --outdir processed/Fairfax/"

#Remove BWA header because it clashes with the Picard MarkDuplicates
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 1000 --ncores 1 --jobname remove_header --command "python ~/software/utils/bam/bamRemoveBwaHeader.py --indir processed/Fairfax/ --outdir processed/Fairfax/"

#Remove duplicates
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 2200 --ncores 1 --jobname remove_duplicates --command "python ~/software/utils/bam/bamRemoveDuplicates.py --indir processed/Fairfax/ --outdir processed/Fairfax/ --execute True --insuffix .reheadered.bam"

#Index bams
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt  | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir processed/Fairfax/ --insuffix .no_duplicates.bam --execute True"

#### Create fragment coverage files ####

#Sort bams by name
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 9000 --ncores 1 --jobname bamsort --command "python ~/software/utils/bam/bamSortName.py --indir processed/Fairfax/ --outdir processed/Fairfax/ --insuffix .no_duplicates.bam"

#Convert BAMS to bed file of fragments
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bamToFragmentBed --command "python ~/software/utils/bam/bamToFragmentBed.py --indir processed/Fairfax/ --outdir processed/Fairfax/ --insuffix .sortedByName.bam --outsuffix .fragments.bed.gz --execute True"

#Convert BED to bigwig
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname bed2bigwig --command "python ~/software/utils/coverage/bed2bigwig.py --indir processed/Fairfax/ --outdir processed/Fairfax/ --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --execute True"

#Clean-up
rm processed/Fairfax/*/*.filtered.bam
rm processed/Fairfax/*/*_2??.bam
rm processed/Fairfax/*/*.sortedByCoord.bam*
rm processed/Fairfax/*/*.reheadered.bam
rm processed/Fairfax/*/*.new_header.txt

#Call narrow and broad peaks
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakNarrow --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Fairfax/ --outdir processed/Fairfax/ --execute True --control processed/OCallaghan/Input_ctrl/Input_ctrl.no_duplicates.bam --broad False"
cat macrophage-gxe-study/data/chromatin/ChIP/Fairfax_samples.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname chipMacsCallPeakBroad --command "python ~/software/utils/coverage/chipMacsCallPeak.py --indir processed/Fairfax/ --outdir processed/Fairfax/ --execute True --control processed/OCallaghan/Input_ctrl/Input_ctrl.no_duplicates.bam --broad True"


