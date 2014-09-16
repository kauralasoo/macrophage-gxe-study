#20/08/14 - Download raw data from IRODS
python ~/software/utils/fetch-irods.py --runid 13569 --laneids 5,6,7 --sampleids 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --dir raw
python ~/software/utils/fetch-irods.py --runid 13518 --laneids 8 --sampleids 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --dir raw
python ~/software/utils/renameFiles.py --samples macrophage-gxe-study/data/DF_rename_runs_patch1.txt --filedir raw --suffix .bam --execute True
bsub -G team170 -o FarmOut/merge_bams.%J.txt "python ~/software/utils/merge-bams.py --runid 13569 --laneids 5,6,7,8 --sampleids 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --indir raw --outdir bams_merged"
python ~/software/utils/renameFiles.py --samples macrophage-gxe-study/data/DF_patch1_samples.txt --suffix .bam --filedir bams_merged/ --execute True

#Sort bams by name
cut -f2 macrophage-gxe-study/data/DF_patch1_samples.txt | python ~/software/utils/submitJobs.py --MEM 20000 --jobname sort_bams --command "python ~/software/utils/sort-bams.py --indir bams_merged/ --outdir bams_nsort/"
cut -f2 macrophage-gxe-study/data/DF_patch1_samples.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname bam2fastq --command "python ~/software/utils/bam2fastq.py --indir bams_nsort/ --outdir fq/"

#Build tophat index of the transcriptome
bsub -G team170 -R "select[mem>6000] rusage[mem=6000]" -M 6000 -o tophat_index.%J.txt "tophat -G genes/Homo_sapiens.GRCh38.76.gtf --transcriptome-index genes/GRCh38_76/GRCh38_76 bowtie2-index/GRCh38"

#Align reads to the reference using Tophat2
cut -f2 macrophage-gxe-study/data/DF_patch1_samples.txt | python ~/software/utils/tophat-align.py --index ../../annotations/GRCh38/bowtie2-index/GRCh38 --txindex ../../annotations/GRCh38/genes/GRCh38_76/GRCh38_76 --out tophat2/ --ncores 6 --fastq fq/ --no_execute True
cut -f2 macrophage-gxe-study/data/DF_patch1_samples.txt | head -n 15 | python ~/software/utils/submitJobs.py --ncores 6 --MEM 8000 --jobname tophat2_align --command "python ~/software/utils/tophat-align.py --index ../../annotations/GRCh38/bowtie2-index/GRCh38 --txindex ../../annotations/GRCh38/genes/GRCh38_76/GRCh38_76 --out tophat2/ --ncores 6 --fastq fq/"

#Some ran out of time. Use the long queue to do the alignment
cut -f1 long.txt | python ~/software/utils/submitJobs.py --ncores 6 --queue long --MEM 8000 --jobname tophat2_align_long --command "python ~/software/utils/tophat-align.py --index ../../annotations/GRCh38/bowtie2-index/GRCh38 --txindex ../../annotations/GRCh38/genes/GRCh38_76/GRCh38_76 --out tophat2/ --ncores 6 --fastq fq/"


