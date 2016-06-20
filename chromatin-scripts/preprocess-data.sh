#Convert Ivashkiv data to fastq
cut -f2 macrophage-chromatin/data/Ivashkiv_2013.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname sra2fastq --command "python ~/software/utils/sra2fastq.py --indir data/Ivashkiv_2013/sra/ --outdir data/Ivashkiv_2013/fastq/"

#Merge fastqs for Rehli_2012 data set
bsub -G team170  -R "span[hosts=1] select[mem>200] rusage[mem=200]" -M 200 -o FarmOut/merge-fastq.%J.jobout "cat macrophage-chromatin/data/Rehli_2012_exp_run.txt | python ~/software/utils/merge-fastq.py --indir data/Rehli_2012/fastq/ --outdir data/Rehli_2012/fastq/"

#Merge fastqs for Rehli_2013 data set
bsub -G team170  -R "span[hosts=1] select[mem>200] rusage[mem=200]" -M 200 -o FarmOut/merge-fastq.%J.jobout "cat macrophage-chromatin/data/Rehli_2013_exp_run.txt | python ~/software/utils/merge-fastq.py --indir data/Rehli_2013/fastq/ --outdir data/Rehli_2013/fastq/"

#Align Rehli_2013 dataset to GRCh38
cut -f2  macrophage-chromatin/data/Rehli_2013.txt |  python ~/software/utils/submitJobs.py --MEM 7000 --ncores 4 --jobname bowtie2-align --command "python ~/software/utils/bowtie2-align.py --fastq data/Rehli_2013/fastq/ --output data/Rehli_2013/bams --index ../../annotations/GRCh38/bowtie2-index/GRCh38 --single_end True --suffix .fastq.gz"

#Align Rehli_2012 dataset to GRCh38 using bowtie2
cut -f2  macrophage-chromatin/data/Rehli_2012.txt |  python ~/software/utils/submitJobs.py --MEM 7000 --ncores 4 --jobname bowtie2-align-rehli2012 --command "python ~/software/utils/bowtie2-align.py --fastq data/Rehli_2012/fastq/ --output data/Rehli_2012/bams --index ../../annotations/GRCh38/bowtie2-index/GRCh38 --single_end True --suffix .fastq.gz"

#Align Rehli_2013 dataset to GRCh38 using bowtie1
cut -f2  macrophage-chromatin/data/Rehli_2013.txt |  python ~/software/utils/submitJobs.py --MEM 4000 --ncores 4 --jobname bowtie1-align-rehli2013 --command "python ~/software/utils/bowtie-align.py --fastq data/Rehli_2013/fastq/ --output data/Rehli_2013/bams_bowtie --index ../../annotations/GRCh38/bowtie-index/GRCh38 --single_end True --suffix .fastq.gz"
