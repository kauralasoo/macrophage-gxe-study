#Fetch all file names from iRODS
python ~/software/utils/irodsGetSamplesInStudy.py --studyName "Genetics of gene expression in human macrophage response to Salmonella" |  cut -f1 -d "." | uniq > fastq/SL1344_samples_4.txt

#Match file names to sample names
python ~/software/utils/irodsFetchMeta.py --irodsList fastq/SL1344_samples.txt | sort -k1 > fastq/SL1344_names.txt 
python ~/software/utils/irodsFetchMeta.py --irodsList fastq/SL1344_samples_2.txt | sort -k1 > fastq/SL1344_names_2.txt 
python ~/software/utils/irodsFetchMeta.py --irodsList fastq/SL1344_samples_3.txt | sort -k1 > fastq/SL1344_names_3.txt 

#Fetch lanelets in cram format from irods
cut -f1 fastq/acLDL_samples.txt | python ~/software/utils/fetch-irods.py --dir fastq/SL1344/cram/ --suffix .cram
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_samples_2.txt | python ~/software/utils/fetch-irods.py --dir fastq/SL1344/cram/ --suffix .cram
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_samples_3.txt | python ~/software/utils/fetch-irods.py --dir fastq/SL1344/cram/ --suffix .cram

#Convert remaining crams to fastq
cut -f1 fastq/remaining.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/bam/cramToBam.py --inputDir fastq/SL1344/ --outputDir fastq/SL1344/"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_samples_2.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/bam/cramToFastq.py --inputDir fastq/SL1344/cram/ --outputDir fastq/SL1344/cram/"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_samples_3.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/bam/cramToFastq.py --inputDir fastq/SL1344/cram/ --outputDir fastq/SL1344/cram/"

#Merge split bams into joint ones and rename
cat fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/ --outdir fastq/SL1344/ --suffix .1.fastq.gz"
cat fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/ --outdir fastq/SL1344/ --suffix .2.fastq.gz"

cat macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/cram/ --outdir fastq/SL1344/ --suffix .1.fastq.gz"
cat macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/cram/ --outdir fastq/SL1344/ --suffix .2.fastq.gz"

cat macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/cram/ --outdir fastq/SL1344/ --suffix .1.fastq.gz"
cat macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/merge-fastq.py --indir fastq/SL1344/cram/ --outdir fastq/SL1344/ --suffix .2.fastq.gz"

#Align reads to the transcriptome using STAR
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR/SL1344/ --fastqDir fastq/SL1344/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_2.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align_160715 --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR/SL1344/ --fastqDir fastq/SL1344/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_3.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align_160715 --ncores 8 --queue hugemem --command "python ~/software/utils/STAR-align.py --outputDir STAR/SL1344/ --fastqDir fastq/SL1344/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"

#Count the number of reads aligning to the mycoplasma genome
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname mycoplasmaTest --command "python ~/software/utils/mycoplasmaTest.py --inputDir fastq/SL1344/ --outdir STAR/SL1344/ --bwaIndex ../../annotations/Mycoplasma/bwa_index/Mycoplasma_genomes.fa --execute True"

#Count reads overlaping GENCODE basic annotations
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"

#Count reads overlapping full Ensembl 79 annotations
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gtf --strand 2 --countsSuffix .counts.txt --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gtf --strand 2 --countsSuffix .counts.txt --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/SL1344/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gtf --strand 2 --countsSuffix .counts.txt --execute True"

#Convert bedgraph files to bigwig
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/bedgraph2bigwig.py --indir STAR/SL1344 --outdir STAR/SL1344 --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

#Index bam files
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --execute True"

#LiftOver raw genotypes from GRCh37 to GRCh38
echo "hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes" | python ~/software/utils/submitJobs.py --MEM 5000 --jobname liftOverVCF --command "python ~/software/utils/vcf/liftoverVcfGenotypes.py --chrMapFwd macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt --chrMapRev macrophage-gxe-study/data/liftOver_genotypes/Hg38ToGRCh38_chromosome_map.txt --liftOver macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain --reference ../../annotations/hg38/hg38.fa --vcfSuffix .vcf.gz --indir genotypes/raw/gtarray/vcf/20141111_858samples/ --outdir genotypes/GRCh38/ --execute True"

#Extract samples from the large VCF file
echo "hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname filterVCF --command "python ~/software/utils/vcf/filterVcf.py  --sampleList macrophage-gxe-study/data/sample_lists/SL1344/SL1344_gt_list.txt --MAF 0.05 --indir genotypes/GRCh38/genotyped/ --outdir genotypes/SL1344/ --execute True  --vcfSuffix .vcf.gz"

#Use verifyBamID to check concordance with the vcf file
cut -f1 fastq/SL1344_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/selected_genotypes.GRCh38.sorted.vcf.gz --execute True" 
cut -f1 fastq/SL1344_failed.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/selected_genotypes.GRCh38.sorted.vcf.gz --execute True" 
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/SL1344/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted.filtered.vcf.gz  --execute True" 
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/SL1344 --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/SL1344/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted.filtered.vcf.gz  --execute True" 

#Count intronic and exonic reads
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | head -n1 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname exonCounts --command "python ~/software/utils/bam/bamToExonCounts.py --sampleDir STAR/SL1344/ --gtf annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.exons.gff3 --strand 2 --countsSuffix .exon_counts.txt --type exon --bamSuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | tail -n 30 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname exonCounts --command "python ~/software/utils/bam/bamToExonCounts.py --sampleDir STAR/SL1344/ --gtf annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.exons.gff3 --strand 2 --countsSuffix .exon_counts.txt --type exon --bamSuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | grep _A | tail -n 20 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname exonCounts --command "python ~/software/utils/bam/bamToExonCounts.py --sampleDir STAR/SL1344/ --gtf annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.exons.gff3 --strand 2 --countsSuffix .exon_counts.txt --type exon --bamSuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | grep _A | tail -n 62 | head -n42 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname exonCounts --command "python ~/software/utils/bam/bamToExonCounts.py --sampleDir STAR/SL1344/ --gtf annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.exons.gff3 --strand 2 --countsSuffix .exon_counts.txt --type exon --bamSuffix .Aligned.sortedByCoord.out.bam --execute True"

cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | head -n1 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname intronCounts --command "python ~/software/utils/bam/bamToExonCounts.py --sampleDir STAR/SL1344/ --gtf annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.introns.gff3 --strand 2 --countsSuffix .intron_counts.txt --type intron --bamSuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | grep _A | tail -n 62 | head -n42 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname intronCounts --command "python ~/software/utils/bam/bamToExonCounts.py --sampleDir STAR/SL1344/ --gtf annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.introns.gff3 --strand 2 --countsSuffix .intron_counts.txt --type intron --bamSuffix .Aligned.sortedByCoord.out.bam --execute True"

#Construct intron events from the intron and exon counts
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | grep _A | python ~/software/utils/bam/constructIntronEvents.py --intronGFF annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.introns.gff3 --exonGFF annotations/exon_intron_annot/Homo_sapiens.GRCh38.79.exons.gff3 --sampleDir STAR/SL1344/

#### PEER ####
#Run PEER on each condition separately using only expressed genes
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python macrophage-gxe-study/SL1344/eQTL/runPEER.py --input results/SL1344/PEER/input/cond_A_exprs.txt  --outdir results/SL1344/PEER/expressed_A_10/ --n_factors 10"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python macrophage-gxe-study/SL1344/eQTL/runPEER.py --input results/SL1344/PEER/input/cond_B_exprs.txt  --outdir results/SL1344/PEER/expressed_B_10/ --n_factors 10"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python macrophage-gxe-study/SL1344/eQTL/runPEER.py --input results/SL1344/PEER/input/cond_C_exprs.txt  --outdir results/SL1344/PEER/expressed_C_10/ --n_factors 10"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python macrophage-gxe-study/SL1344/eQTL/runPEER.py --input results/SL1344/PEER/input/cond_D_exprs.txt  --outdir results/SL1344/PEER/expressed_D_10/ --n_factors 10"

#Run PEER on each condition separately using all genes
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python macrophage-gxe-study/SL1344/eQTL/runPEER.py --input results/SL1344/PEER/input/cond_A_exprs.all_genes.txt  --outdir results/SL1344/PEER/all_A_10/ --n_factors 10"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python macrophage-gxe-study/SL1344/eQTL/runPEER.py --input results/SL1344/PEER/input/cond_B_exprs.all_genes.txt  --outdir results/SL1344/PEER/all_B_10/ --n_factors 10"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python macrophage-gxe-study/SL1344/eQTL/runPEER.py --input results/SL1344/PEER/input/cond_C_exprs.all_genes.txt  --outdir results/SL1344/PEER/all_C_10/ --n_factors 10"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python macrophage-gxe-study/SL1344/eQTL/runPEER.py --input results/SL1344/PEER/input/cond_D_exprs.all_genes.txt  --outdir results/SL1344/PEER/all_D_10/ --n_factors 10"





