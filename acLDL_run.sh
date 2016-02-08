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
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python ~/software/utils/runPEER.py --input results/acLDL/PEER/input//Ctrl.exprs.txt --outdir results/acLDL/PEER/Ctrl_10/ --n_factors 10"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>500] rusage[mem=500]" -q normal -M 500 -o FarmOut/PEER.%J.jobout "python ~/software/utils/runPEER.py --input results/acLDL/PEER/input//AcLDL.exprs.txt --outdir results/acLDL/PEER/AcLDL_10/ --n_factors 10"


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


