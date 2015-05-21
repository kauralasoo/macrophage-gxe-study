#Fetch list of lanelets from IRODS
python ~/software/utils/irodsGetSamplesInStudy.py --studyName "Characterization of iPSC derived macrophages - cardiovascular pilot" |  cut -f1 -d "." | uniq > fastq/acLDL_samples.txt

#Map lanelet ids to sample ids
python ~/software/utils/irodsFetchMeta.py --irodsList fastq/acLDL_samples.txt | sort -k1 > fastq/acLDL_names.txt 

#Fetch lanelets in cram format from irods
cut -f1 fastq/acLDL_samples.txt | python ~/software/utils/fetch-irods.py --dir fastq/acLDL/cram/ --suffix .cram

#Convert cram files into fastq
cut -f1 fastq/acLDL_samples.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/cramToFastq.py --inputDir fastq/acLDL/cram/ --outputDir fastq/acLDL/"

#Merge split fastq into joint ones and rename
cat fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/acLDL/ --outdir fastq/acLDL/ --suffix .1.fastq.gz"
cat fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/acLDL/ --outdir fastq/acLDL/ --suffix .2.fastq.gz"

#Align reads to the genome using STAR
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR/acLDL/ --fastqDir fastq/acLDL/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"

#Count reads overlaping GENCODE basic annotations
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/acLDL/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"

#Count reads overlapping full Ensembl 79 annotations
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/acLDL/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gtf --strand 2 --countsSuffix .counts.txt --execute True"

#Convert bedgraph to bigwig and compress bedgraphs
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"