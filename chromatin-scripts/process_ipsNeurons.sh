#Convert .cram to BAM
samtools view -b BULB_1_ATAC/BULB_1_ATAC.cram > BULB_1_ATAC/BULB_1_ATAC.bam
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 -o bwa_index.%J.jobout "samtools view -b CUPS_3_ATAC/CUPS_3_ATAC.cram > CUPS_3_ATAC/CUPS_3_ATAC.bam"

#Sort bams by name
cut -f1 data/ipsNeurons/sample_names.txt | python ~/software/utils/submitJobs.py --MEM 21000 --jobname sort_bams --command "python ~/software/utils/sort-bams.py --indir data/ipsNeurons --outdir data/ipsNeurons --insuffix .bam" 

#Convert BAM to fragment BED
cut -f1 data/ipsNeurons/sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname sort_bams --command "python ~/software/utils/bamToFragmentBed.py --indir data/ipsNeurons --outdir data/ipsNeurons --insuffix .sortedByName.bam --outsuffix .fragments.bed.gz --execute True"

#Convert BED to BigWig
cut -f1 data/ipsNeurons/sample_names.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname sort_bams --command "python ~/software/utils/bed2bigwig.py --indir data/ipsNeurons --outdir data/ipsNeurons --chrlengths data/ipsNeurons/header.txt --execute True"
