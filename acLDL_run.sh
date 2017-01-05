#Fetch list of lanelets from IRODS
python ~/software/utils/irods/irodsGetSamplesInStudy.py --studyName "Characterization of iPSC derived macrophages - cardiovascular pilot" |  cut -f1 -d "." | uniq > fastq/acLDL_samples.txt

#Fetch and updated list of lanelets
python ~/software/utils/irods/irodsGetSamplesInStudy.py --studyName "Characterization of iPSC derived macrophages - cardiovascular pilot" |  cut -f1 -d "." | uniq > fastq/acLDL_samples_2.txt
python ~/software/utils/irods/irodsGetSamplesInStudy.py --studyName "Characterization of iPSC derived macrophages - cardiovascular pilot" |  cut -f1 -d "." | uniq > fastq/acLDL_samples_3.txt

#Map lanelet ids to sample ids
python ~/software/utils/irods/irodsFetchMeta.py --irodsList fastq/acLDL_samples.txt | sort -k1 > fastq/acLDL_names.txt 
python ~/software/utils/irods/irodsFetchMeta.py --irodsList fastq/acLDL_samples_2.txt | sort -k1 > fastq/acLDL_names_2.txt 
python ~/software/utils/irods/irodsFetchMeta.py --irodsList fastq/acLDL_samples_3.txt | sort -k1 > fastq/acLDL_names_3.txt 
python ~/software/utils/irods/irodsFetchMeta.py --irodsList fastq/acLDL_samples_3.txt | sort -k1 > fastq/acLDL_names_4.txt 

#Fetch lanelets in cram format from irods
cut -f1 fastq/acLDL_samples.txt | python ~/software/utils/irods/fetch-irods.py --dir fastq/acLDL/cram/ --suffix .cram
cut -f1 fastq/acLDL_samples_2.txt | python ~/software/utils/irods/fetch-irods.py --dir fastq/acLDL/cram/ --suffix .cram
cut -f1 fastq/acLDL_samples_3.txt | python ~/software/utils/irods/fetch-irods.py --dir fastq/acLDL/cram/ --suffix .cram
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_samples_5.txt | python ~/software/utils/irods/fetch-irods.py --dir fastq/acLDL/cram/ --suffix .cram

#Convert cram files into fastq
cut -f1 fastq/acLDL_samples.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/cramToFastq.py --inputDir fastq/acLDL/cram/ --outputDir fastq/acLDL/"
cut -f1 fastq/acLDL_samples_2.txt |  python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/cramToFastq.py --inputDir fastq/acLDL/cram/ --outputDir fastq/acLDL/"
cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_samples_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/bam/cramToFastq.py --inputDir fastq/acLDL/cram/ --outputDir fastq/acLDL/cram/"
cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_samples_4.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/bam/cramToFastq.py --inputDir fastq/acLDL/cram/ --outputDir fastq/acLDL/cram/"
cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_samples_5.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname cramToFastq --command "python ~/software/utils/bam/cramToFastq.py --inputDir fastq/acLDL/cram/ --outputDir fastq/acLDL/cram/"

#Merge split fastq into joint ones and rename
cat fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/acLDL/ --outdir fastq/acLDL/ --suffix .1.fastq.gz"
cat fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/acLDL/ --outdir fastq/acLDL/ --suffix .2.fastq.gz"

cat fastq/acLDL_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/acLDL/ --outdir fastq/acLDL/ --suffix .1.fastq.gz"
cat fastq/acLDL_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/acLDL/ --outdir fastq/acLDL/ --suffix .2.fastq.gz"

cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/acLDL/cram/ --outdir fastq/acLDL/ --suffix .1.fastq.gz"
cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/acLDL/cram/ --outdir fastq/acLDL/ --suffix .2.fastq.gz"

cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_4.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/acLDL/cram/ --outdir fastq/acLDL/ --suffix .1.fastq.gz"
cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_4.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/acLDL/cram/ --outdir fastq/acLDL/ --suffix .2.fastq.gz"

cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_5.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/acLDL/cram/ --outdir fastq/acLDL/ --suffix .1.fastq.gz"
cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_5.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname merge_bams --command "python ~/software/utils/fastq/merge-fastq.py --indir fastq/acLDL/cram/ --outdir fastq/acLDL/ --suffix .2.fastq.gz"

#Align reads to the genome using STAR
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/align/STAR-align.py --outputDir STAR/acLDL/ --fastqDir fastq/acLDL/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"
cut -f1 fastq/acLDL_names_2.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/align/STAR-align.py --outputDir STAR/acLDL/ --fastqDir fastq/acLDL/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_3.txt | python ~/software/utils/submitJobs.py --MEM 32000 --jobname star_align --ncores 8 --queue hugemem --command "python ~/software/utils/align/STAR-align.py --outputDir STAR/acLDL/ --fastqDir fastq/acLDL/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_4.txt | python ~/software/utils/submitJobs.py --MEM 35000 --jobname star_align_051115 --ncores 8 --queue hugemem --command "python ~/software/utils/align/STAR-align.py --outputDir STAR/acLDL/ --fastqDir fastq/acLDL/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_5.txt | python ~/software/utils/submitJobs.py --MEM 35000 --jobname star_align_230116 --ncores 8 --queue hugemem --command "python ~/software/utils/align/STAR-align.py --outputDir STAR/acLDL/ --fastqDir fastq/acLDL/ --genomeDir ../../annotations/GRCh38/STAR_index_79/ --runThreadN 8"


#Count reads overlaping GENCODE basic annotations
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir STAR/acLDL/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"
cut -f1 fastq/acLDL_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir STAR/acLDL/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir STAR/acLDL/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_4.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir STAR/acLDL/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_5.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam/bam2counts.py --sampleDir STAR/acLDL/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gencode_basic.gtf --strand 2 --countsSuffix .gencode_basic.counts.txt --execute True --donotsort True"

#Count reads overlapping full Ensembl 79 annotations
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/acLDL/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gtf --strand 2 --countsSuffix .counts.txt --execute True"
cut -f1 fastq/acLDL_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/acLDL/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gtf --strand 2 --countsSuffix .counts.txt --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname featureCounts --command "python ~/software/utils/bam2counts.py --sampleDir STAR/acLDL/ --gtf ../../annotations/GRCh38/genes/Homo_sapiens.GRCh38.79.gtf --strand 2 --countsSuffix .counts.txt --execute True"

#Convert bedgraph to bigwig and compress bedgraphs
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/coverage/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/coverage/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

cut -f1 fastq/acLDL_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/coverage/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 fastq/acLDL_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/coverage/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/coverage/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/coverage/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_4.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/coverage/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_4.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/coverage/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_5.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/coverage/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str1.out.bg --outsuffix .str1.bw"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_5.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedgraph2bigwig --command "python ~/software/utils/coverage/bedgraph2bigwig.py --indir STAR/acLDL --outdir STAR/acLDL --chrlengths ../../annotations/GRCh38/bt2-index/chromosome_lengths.txt --insuffix .Signal.Unique.str2.out.bg --outsuffix .str2.bw"

#Construct genotype list
cut -f2 acLDL_sample_genotype_map.txt | sort | uniq > acLDL_genotype_list.txt

#Extract samples from the large VCF file
echo "hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname filterVCF --command "python ~/software/utils/vcf/filterVcf.py  --sampleList macrophage-gxe-study/data/sample_lists/acLDL/acLDL_gt_list.txt --MAF 0.05 --indir genotypes/GRCh38/genotyped/ --outdir genotypes/acLDL/ --execute True  --vcfSuffix .vcf.gz"

#Index bam files
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 fastq/acLDL_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_3.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_4.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --execute True"
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_5.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --execute True"

#Run VerifyBamID
cut -f1 fastq/acLDL_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/acLDL/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted.filtered.vcf.gz --execute True" 
cut -f1 fastq/acLDL_names_2.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/acLDL/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted.filtered.vcf.gz --execute True" 
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_all.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/acLDL/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted.filtered.vcf.gz --execute True" 
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_4.txt | python ~/software/utils/submitJobs.py --MEM 1500 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/acLDL/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted.filtered.vcf.gz --execute True" 
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_names_5.txt | python ~/software/utils/submitJobs.py --MEM 1500 --jobname verifyBamID --command  "python ~/software/utils/runVerifyBamID.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --vcf genotypes/acLDL/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.GRCh38.sorted.filtered.vcf.gz --execute True" 

#### PEER ####
#Run PEER on each condition separately using only expressed genes
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python ~/software/utils/runPEER.py --input results/acLDL/PEER/input//Ctrl.exprs.txt --outdir results/acLDL/PEER/output/Ctrl --n_factors 10"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python ~/software/utils/runPEER.py --input results/acLDL/PEER/input//AcLDL.exprs.txt --outdir results/acLDL/PEER/output/AcLDL --n_factors 10"


#### Count ASE ####
#Add read group to the BAM files to make them work with GATK
cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_sample_gt_map.txt | head -n 100 | tail -n 99 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamAddRG --command "python ~/software/utils/bam/bamAddRG.py --indir STAR/acLDL --outdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --outsuffix .Aligned.sortedByCoord.RG.bam --execute True"
cat macrophage-gxe-study/data/sample_lists/acLDL/acLDL_sample_gt_map.txt | tail -n 50 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bamAddRG --command "python ~/software/utils/bam/bamAddRG.py --indir STAR/acLDL --outdir STAR/acLDL --insuffix .Aligned.sortedByCoord.out.bam --outsuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Index bams
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_sample_gt_map.txt | tail -n 149 | python ~/software/utils/submitJobs.py --MEM 1000 --jobname index_bams --command  "python ~/software/utils/bam/index-bams.py --bamdir STAR/acLDL --insuffix .Aligned.sortedByCoord.RG.bam --execute True"

#Index the vcf file
tabix -p vcf genotypes/acLDL/imputed_20151005/imputed.70_samples.snps_only.vcf.gz
bcftools index genotypes/acLDL/imputed_20151005/imputed.70_samples.snps_only.vcf.gz

#Use ASEReadCounter to count allele-specific expression
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_sample_gt_map.txt | python ~/software/utils/submitJobs.py --MEM 10000 --jobname bamCountASE --command "python ~/software/utils/rasqual/bamCountASE.py --indir STAR/acLDL --outdir STAR/acLDL --insuffix .Aligned.sortedByCoord.RG.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/acLDL/imputed_20151005/imputed.70_samples.snps_only.vcf.gz --execute True --Xmx 8g"

#Construct a sample-sample map for meregeASECounts.py script
cut -f1 macrophage-gxe-study/data/sample_lists/acLDL/acLDL_sample_gt_map.txt | awk -v OFS='\t' '{print $1, $1}' > results/acLDL/rasqual/input/sample_sample_map.txt

#Merge all allele-specific counts into one matrix
echo "mergeASECounts" | python ~/software/utils/submitJobs.py --MEM 32000 --jobname mergeASECounts --queue hugemem --command "python ~/software/utils/rasqual/mergeASECounts.py --sample_list results/acLDL/rasqual/input/sample_sample_map.txt --indir STAR/acLDL --suffix .ASEcounts > results/acLDL/combined_ASE_counts.txt"

#Sort and index the ASE counts file
(zcat results/acLDL/combined_ASE_counts.txt.gz | head -n1 && zcat results/acLDL/combined_ASE_counts.txt.gz | tail -n+2 | sort -k1,1 -k2,2n) | bgzip > results/acLDL/combined_ASE_counts.sorted.txt.gz
tabix -s 1 -b 2 -e 2 -S 1 results/acLDL/combined_ASE_counts.sorted.txt.gz 


#Extract genotype ids for each condition
cut -f2 results/acLDL/rasqual/input/Ctrl.sg_map.txt > results/acLDL/rasqual/input/Ctrl.genotypes.txt
cut -f2 results/acLDL/rasqual/input/AcLDL.sg_map.txt > results/acLDL/rasqual/input/AcLDL.genotypes.txt

#Extract samples from the global VCF file
bcftools view -Oz -S results/acLDL/rasqual/input/AcLDL.genotypes.txt genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/acLDL/rasqual/input/AcLDL.vcf &
bcftools view -Oz -S results/acLDL/rasqual/input/Ctrl.genotypes.txt genotypes/acLDL/imputed_20151005/imputed.70_samples.sorted.filtered.named.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - > results/acLDL/rasqual/input/Ctrl.vcf &

#Add ASE counts into the VCF file
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 32000 --jobname vcfAddASE --queue hugemem --command "python ~/software/rasqual/scripts/vcfAddASE.py --ASEcounts results/acLDL/combined_ASE_counts.txt --ASESampleGenotypeMap results/acLDL/rasqual/input/Ctrl.sg_map.txt --VCFfile results/acLDL/rasqual/input/Ctrl.vcf | bgzip > results/acLDL/rasqual/input/Ctrl.ASE.vcf.gz"
echo "vcfAddASE" | python ~/software/utils/submitJobs.py --MEM 32000 --jobname vcfAddASE --queue hugemem --command "python ~/software/rasqual/scripts/vcfAddASE.py --ASEcounts results/acLDL/combined_ASE_counts.txt --ASESampleGenotypeMap results/acLDL/rasqual/input/AcLDL.sg_map.txt --VCFfile results/acLDL/rasqual/input/AcLDL.vcf | bgzip > results/acLDL/rasqual/input/AcLDL.ASE.vcf.gz"

#Index VCF files
tabix -p vcf results/acLDL/rasqual/input/Ctrl.ASE.vcf.gz 
tabix -p vcf results/acLDL/rasqual/input/AcLDL.ASE.vcf.gz 

#RASQUAL (no covariates, library sizes + GC), 100kb
cat results/acLDL/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/acLDL/rasqual/input/Ctrl.expression.bin  --offsets results/acLDL/rasqual/input/Ctrl.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/Ctrl.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_100kb.txt --outprefix results/acLDL/rasqual/output/chr11_naive_100kb_gc --execute True"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/utils/rasqual/mergeRasqualBatches.py --prefix results/acLDL/rasqual/output/chr11_naive_100kb_gc"

#RASQUAL (SVD covariates, library sizes + GC), 100kb
cat results/acLDL/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/acLDL/rasqual/input/Ctrl.expression.bin  --offsets results/acLDL/rasqual/input/Ctrl.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/Ctrl.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_100kb.txt --outprefix results/acLDL/rasqual/output/chr11_naive_100kb_gc_svd --covariates results/acLDL/rasqual/input/Ctrl.svd_covariates.bin --execute True"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/utils/rasqual/mergeRasqualBatches.py --prefix results/acLDL/rasqual/output/chr11_naive_100kb_gc_svd"

#RASQUAL (4 PEER covariates, library sizes + GC), 100kb
cat results/acLDL/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/acLDL/rasqual/input/Ctrl.expression.bin  --offsets results/acLDL/rasqual/input/Ctrl.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/Ctrl.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_100kb.txt --outprefix results/acLDL/rasqual/output/chr11_naive_100kb_gc_PEER --covariates results/acLDL/rasqual/input/Ctrl.PEER_covariates.bin --execute True"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/utils/rasqual/mergeRasqualBatches.py --prefix results/acLDL/rasqual/output/chr11_naive_100kb_gc_PEER"

#RASQUAL (2 PEER covariates + sex, library sizes + GC), 500kb
cat results/acLDL/rasqual/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_Ctrl --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/acLDL/rasqual/input/Ctrl.expression.bin  --offsets results/acLDL/rasqual/input/Ctrl.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/Ctrl.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_500kb.txt --outprefix results/acLDL/rasqual/output/Ctrl_500kb/batches/Ctrl_500kb --covariates results/acLDL/rasqual/input/Ctrl.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force' --execute True"

cat results/acLDL/rasqual/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname runRasqual_AcLDL --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/acLDL/rasqual/input/AcLDL.expression.bin  --offsets results/acLDL/rasqual/input/AcLDL.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/AcLDL.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_500kb.txt --outprefix results/acLDL/rasqual/output/AcLDL_500kb/batches/AcLDL_500kb --covariates results/acLDL/rasqual/input/AcLDL.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force' --execute True"

#Run RASQUAL with a single permutation to calculate empirical FDR
cat results/acLDL/rasqual/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500  --ncores 8 --jobname runRasqual_Ctrl --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/acLDL/rasqual/input/Ctrl.expression.bin  --offsets results/acLDL/rasqual/input/Ctrl.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/Ctrl.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_500kb.txt --outprefix results/acLDL/rasqual/output/Ctrl_500kb/batches/Ctrl_500kb --covariates results/acLDL/rasqual/input/Ctrl.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force --n-threads 8 --random-permutation' --execute True"

cat results/acLDL/rasqual/input/gene_batches.txt | python ~/software/utils/submitJobs.py --MEM 500  --ncores 8 --jobname runRasqual_AcLDL --command "python ~/software/rasqual/scripts/runRasqual.py --readCounts results/acLDL/rasqual/input/AcLDL.expression.bin  --offsets results/acLDL/rasqual/input/AcLDL.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/AcLDL.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_500kb.txt --outprefix results/acLDL/rasqual/output/AcLDL_500kb/batches/AcLDL_500kb --covariates results/acLDL/rasqual/input/AcLDL.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force --n-threads 8 --random-permutation'--execute True"


#Rerun failed samples
cat results/acLDL/rasqual/input/ctrl_failed_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --queue normal --jobname runRasqual_Ctrl3 --ncores 4 --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/acLDL/rasqual/input/Ctrl.expression.bin  --offsets results/acLDL/rasqual/input/Ctrl.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/Ctrl.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_500kb.txt --outprefix results/acLDL/rasqual/output/Ctrl_500kb/batches/Ctrl_500kb --covariates results/acLDL/rasqual/input/Ctrl.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force --n-threads 4' --execute True"

cat results/acLDL/rasqual/input/acldl_failed_batches.txt | python ~/software/utils/submitJobs.py --MEM 500 --queue normal --jobname runRasqual_AcLDL3 --ncores 4 --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/acLDL/rasqual/input/AcLDL.expression.bin  --offsets results/acLDL/rasqual/input/AcLDL.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/AcLDL.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_500kb.txt --outprefix results/acLDL/rasqual/output/AcLDL_500kb/batches/AcLDL_500kb --covariates results/acLDL/rasqual/input/AcLDL.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force --n-threads 4' --execute True"

#Rerun the remanining failed samples using 10 threads and stricter IMP2 quality score filtering
cat ctrl_failed_batches_2.txt | python ~/software/utils/submitJobs.py --MEM 500 --queue normal --jobname runRasqual_Ctrl2 --ncores 5 --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/acLDL/rasqual/input/Ctrl.expression.bin  --offsets results/acLDL/rasqual/input/Ctrl.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/Ctrl.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_500kb.txt --outprefix results/acLDL/rasqual/output/Ctrl_500kb/batches/Ctrl_500kb --covariates results/acLDL/rasqual/input/Ctrl.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force --n-threads 5 --imputation-quality-fsnp 0.9' --execute True"

cat acLDL_failed_batches_2.txt | python ~/software/utils/submitJobs.py --MEM 500 --queue long --jobname runRasqual_AcLDL2 --ncores 2 --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/acLDL/rasqual/input/AcLDL.expression.bin  --offsets results/acLDL/rasqual/input/AcLDL.gc_library_size.bin --n 70 --geneids results/acLDL/rasqual/input/feature_names.txt --vcf results/acLDL/rasqual/input/AcLDL.ASE.vcf.gz --geneMetadata results/acLDL/rasqual/input/gene_snp_count_500kb.txt --outprefix results/acLDL/rasqual/output/AcLDL_500kb/batches/AcLDL_500kb --covariates results/acLDL/rasqual/input/AcLDL.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--force --n-threads 2 --imputation-quality-fsnp 0.9' --execute True"

#Merge batches
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/acLDL/rasqual/output/Ctrl_500kb/batches/Ctrl_500kb"
echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/rasqual/scripts/mergeRasqualBatches.py --prefix results/acLDL/rasqual/output/AcLDL_500kb/batches/AcLDL_500kb"

#Move merged file to parent dir
mv results/acLDL/rasqual/output/Ctrl_500kb/batches/Ctrl_500kb.txt results/acLDL/rasqual/output/Ctrl_500kb/
mv results/acLDL/rasqual/output/AcLDL_500kb/batches/AcLDL_500kb.txt results/acLDL/rasqual/output/AcLDL_500kb/

# Sort and filter merged p-values
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.txt | sort -k3,3 -k4,4n | bgzip > results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.sorted.txt.gz"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/sortRasqual.%J.jobout "grep -v SKIPPED results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.txt | sort -k3,3 -k4,4n | bgzip > results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.sorted.txt.gz"

#Index the output files using Tabix
tabix -s3 -b4 -e4 -f results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.sorted.txt.gz
tabix -s3 -b4 -e4 -f results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.sorted.txt.gz

#Identify completed genes
cut -f1 results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.txt | uniq > results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.completed_ids.txt
cut -f1 results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.txt | uniq > results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.completed_ids.txt

#Compress batches 
tar czf results/acLDL/rasqual/output/AcLDL_500kb/batches.tar.gz results/acLDL/rasqual/output/AcLDL_500kb/batches/ &
tar czf results/acLDL/rasqual/output/Ctrl_500kb/batches.tar.gz results/acLDL/rasqual/output/Ctrl_500kb/batches/ &

### eigenMT ####
#Convert rasqual output into format suitable for eigenMT
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.txt > results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.eigenMT_input.txt"
echo "rasqualToEigenMT" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname rasqualToEigenMT --command "python ~/software/rasqual/scripts/rasqualToEigenMT.py --rasqualOut results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.txt > results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.eigenMT_input.txt"

#Split vcf into chromosomes
cat ../../../macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/vcf/vcfSplitByChromosome.py --vcf imputed.86_samples.sorted.filtered.named.INFO_07.vcf.gz --outdir chromosomes_INFO_07/ --execute False

#Convert vcfs to GDS
/software/R-3.1.2/bin/Rscript ~/software/utils/vcf/vcfToGds.R --vcf-directory chromosomes_INFO_07 --chr-list ../../../macrophage-gxe-study/data/sample_lists/chromosome_list.txt

#Run eigenMT on on each chromosome
cat macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname eigenMTbyChromosome2   --command "python ~/software/utils/eigenMTbyChromosome.py --genepos results/acLDL/rasqual/input/gene_positions.txt --chromosome_dir genotypes/acLDL/imputed_20151005/chromosomes_INFO_07/ --QTL results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.eigenMT_input.txt --out_prefix results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb --cis_dist 1e7 --eigenMT_path ~/software/eigenMT/eigenMT.py"
cat macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/submitJobs.py --MEM 2000 --jobname eigenMTbyChromosome --command "python ~/software/utils/eigenMTbyChromosome.py --genepos results/acLDL/rasqual/input/gene_positions.txt --chromosome_dir genotypes/acLDL/imputed_20151005/chromosomes_INFO_07/ --QTL results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.eigenMT_input.txt --out_prefix results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb --cis_dist 1e7 --eigenMT_path ~/software/eigenMT/eigenMT.py"

# Concat eigenMT results
cat results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.chr_*.eigenMT.txt | grep -v snps > results/acLDL/rasqual/output/Ctrl_500kb/Ctrl_500kb.eigenMT.txt
cat results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.chr_*.eigenMT.txt | grep -v snps > results/acLDL/rasqual/output/AcLDL_500kb/AcLDL_500kb.eigenMT.txt


