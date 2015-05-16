#Construct STAR index from the fasta and GTF file
bsub -G team170 -n4 -R "span[hosts=1] select[mem>40000] rusage[mem=40000]" -q hugemem -M 40000 -o star_index.%J.jobout "STAR --runThreadN 4 --runMode genomeGenerate --genomeDir STAR_index/ --genomeFastaFiles dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile genes/Homo_sapiens.GRCh38.78.gtf --sjdbOverhang 74"

#Download raw data from irods
cut -f1 fastq/acLDL_samples.txt |  python ~/software/utils/irodsToFastq.py --output fastq/acLDL/

#Rename Run ids to sample names
python ~/software/utils/renameFiles.py --samples fastq/acLDL_samples.txt --filedir fastq/acLDL/ --suffix .1.fastq.gz --execute True
python ~/software/utils/renameFiles.py --samples fastq/acLDL_samples.txt --filedir fastq/acLDL/ --suffix .2.fastq.gz --execute True

#Align reads with STAR
cut -f2 fastq/acLDL_samples.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR --fastqDir fastq/acLDL/ --genomeDir ../../annotations/GRCh38/STAR_index/ --runThreadN 8"

#Count reads over exons
cut -f2 fastq/acLDL_samples.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.78.gtf --strand 2 --execute True"

#Construct STAR index (Ensembl 79) from the fasta and GTF file
bsub -G team170 -n4 -R "span[hosts=1] select[mem>40000] rusage[mem=40000]" -q hugemem -M 40000 -o star_index.%J.jobout "STAR --runThreadN 4 --runMode genomeGenerate --genomeDir STAR_index_79/ --genomeFastaFiles dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile genes/Homo_sapiens.GRCh38.79.gtf --sjdbOverhang 74"
