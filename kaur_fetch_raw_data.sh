#Fetch all file names from iRODS
python ~/software/utils/irodsGetSamplesInStudy.py --studyName "Genetics of gene expression in human macrophage response to Salmonella" |  cut -f1 -d "." | uniq > fastq/SL1344_samples.txt

#Match file names to sample names
python ~/software/utils/irodsFetchMeta.py --irodsList fastq/SL1344_samples.txt | sort -k1 > fastq/SL1344_names.txt 

#Download fastq files from iRODS
cut -f1 fastq/SL1344_samples.txt | head -n 100 |  python ~/software/utils/irodsToFastq.py --output fastq/SL1344/
cut -f1 fastq/SL1344_samples.txt | head -n 200 | tail -n 100 |  python ~/software/utils/irodsToFastq.py --output fastq/SL1344/
cut -f1 fastq/SL1344_samples.txt | tail -n 328 |  python ~/software/utils/irodsToFastq.py --output fastq/SL1344/

#Convert remaining crams to fastq
cut -f1 fastq/remaining.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/cramToBam.py --inputDir fastq/SL1344/ --outputDir fastq/SL1344/"

#Merge split bams into joint ones and rename
cat fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/ --outdir fastq/SL1344/ --suffix .1.fastq.gz"
cat fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/ --outdir fastq/SL1344/ --suffix .2.fastq.gz"

#Align reads to the transcriptome using STAR
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR/SL1344/ --fastqDir fastq/SL1344/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"

#Count the number of reads overlapping gene annotations
cut -f1 fastq/SL1344_names.txt | tail -n 40 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.78.gtf --strand 2 --execute True"

cut -f1 fastq/SL1344_names.txt | head -n 92 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.78.gtf --strand 2 --execute True"

cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname mycoplasmaTest --command "python ~/software/utils/mycoplasmaTest.py --inputDir fastq/SL1344/ --outdir STAR/SL1344/ --bwaIndex ../../annotations/Mycoplasma/bwa_index/Mycoplasma_genomes.fa --execute True"

cut -f1 fastq/SL1344_names.txt | head -n 1 | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR1 --fastqDir fastq/SL1344/ --genomeDir ../../annotations/GRCh38/STAR_index/ --runThreadN 8"
