#Convert cram files into fastq
cat sample_list.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/bam/cramToFastq.py --inputDir cram/ --outputDir fastq/"

#MErge Fastq
cat sample_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/ --outdir fastq/ --suffix .2.fastq.gz"
cat sample_names | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/ --outdir fastq/ --suffix .1.fastq.gz"

#Run Snakemake
snakemake --cluster ../../Blood_ATAC/scripts/snakemake_submit.py -np -s acLDL_alternative_transcription.snakefile processed/macroMap/out.txt --jobs 100 --configfile macroMap/macroMap_config.yaml
