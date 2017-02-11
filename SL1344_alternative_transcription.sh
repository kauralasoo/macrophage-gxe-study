#Construct transcriptome index for Salmon
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o construct_index.%J.jobout "salmon --no-version-check index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i salmon_index_85"

#Combined index of cDNA and ncRNA
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o construct_index.%J.jobout "salmon --no-version-check index -t Homo_sapiens.GRCh38.cdna_ncrna.v87.fa.gz -i salmon_index_87"


#Quantify Ensembl transcript expression using Salmon
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/submitJobs.py --MEM 10000 --jobname salmon_quant --ncores 8 --queue normal --command "python ~/software/utils/align/salmonQuant.py --outputDir processed/SL1344_salmon/ --outputSubdir ensembl_full_bootstraps_50 --fastqDir fastq/SL1344 --index ../../annotations/GRCh38/genes/Ensembl_85/salmon_index_85/ --libType ISR --geneMap ../../annotations/GRCh38/genes/Ensembl_85/Homo_sapiens.GRCh38.85.gtf --nCores 8 --numBootstraps 50"

#Copy gene and transcript estimates back to the main results folder
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/copySalmonOutput.py --currentDir processed/SL1344_salmon/ --currentSubdir ensembl_full_bootstraps_50 --newDir STAR/SL1344/ --filename quant.sf --suffix ensembl85
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/copySalmonOutput.py --currentDir processed/SL1344_salmon/ --currentSubdir ensembl_full_bootstraps_50 --newDir STAR/SL1344/ --filename quant.genes.sf --suffix ensembl85

#Estimate dispersions with DRIMSeq
echo "test" | python ~/software/utils/submitJobs.py --MEM 12000 --jobname drimseq_dispersions --ncores 10 --queue normal --command "/software/R-3.3.0/bin/Rscript macrophage-gxe-study/SL1344/splicing/differentialTranscriptExpression.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname drimseq_dispersions --ncores 1 --queue normal --command "/software/R-3.3.0/bin/Rscript macrophage-gxe-study/SL1344/splicing/differentialTranscriptExpression.R"

#Construct new alternative transcription events
cat annotation_batches.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname reviseAnnotations4 --ncores 1 --queue normal --command "/software/R-3.2.2/bin/Rscript macrophage-gxe-study/munge/constructTranscriptionEvents.R"

#Merge all GFF files together
cat reviseAnnotations_batch_* | grep -v "^#" > reviseAnnotations.gff3

#Convert GFF files to fastq
/software/team82/cufflinks/2.2.1/bin/gffread -w reviseAnnotations.contained.fa -g ../../../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa reviseAnnotations.contained.gff3 
/software/team82/cufflinks/2.2.1/bin/gffread -w reviseAnnotations.upstream.fa -g ../../../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa reviseAnnotations.upstream.gff3
/software/team82/cufflinks/2.2.1/bin/gffread -w reviseAnnotations.downstream.fa -g ../../../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa reviseAnnotations.downstream.gff3

#Construct Salmon indices for each type of events
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o construct_index.%J.jobout "salmon --no-version-check index -t reviseAnnotations.contained.fa -i reviseAnnotations.contained.index"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o construct_index.%J.jobout "salmon --no-version-check index -t reviseAnnotations.downstream.fa -i reviseAnnotations.downstream.index"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o construct_index.%J.jobout "salmon --no-version-check index -t reviseAnnotations.upstream.fa -i reviseAnnotations.upstream.index"

#Run quantification on each index
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/submitJobs.py --MEM 4000 --jobname salmon_quant --ncores 8 --queue normal --command "python ~/software/utils/align/salmonQuant.py --outputDir processed/SL1344_salmon/ --outputSubdir contained --fastqDir fastq/SL1344 --index results/reviseAnnotations/reviseAnnotations.contained.index/ --libType ISR --nCores 8 --numBootstraps 0"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/submitJobs.py --MEM 4000 --jobname salmon_quant --ncores 8 --queue normal --command "python ~/software/utils/align/salmonQuant.py --outputDir processed/SL1344_salmon/ --outputSubdir upstream --fastqDir fastq/SL1344 --index results/reviseAnnotations/reviseAnnotations.upstream.index/ --libType ISR --nCores 8 --numBootstraps 0"
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/submitJobs.py --MEM 4000 --jobname salmon_quant --ncores 8 --queue normal --command "python ~/software/utils/align/salmonQuant.py --outputDir processed/SL1344_salmon/ --outputSubdir downstream --fastqDir fastq/SL1344 --index results/reviseAnnotations/reviseAnnotations.downstream.index/ --libType ISR --nCores 8 --numBootstraps 0"

#Copy expression estimates back to main results folder
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/copySalmonOutput.py --currentDir processed/SL1344_salmon/ --currentSubdir contained --newDir STAR/SL1344/ --filename quant.sf --suffix contained &
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/copySalmonOutput.py --currentDir processed/SL1344_salmon/ --currentSubdir upstream --newDir STAR/SL1344/ --filename quant.sf --suffix upstream &
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/copySalmonOutput.py --currentDir processed/SL1344_salmon/ --currentSubdir downstream --newDir STAR/SL1344/ --filename quant.sf --suffix downstream &

#Merge all event types into a signle file
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/mergeSalmonEstimates.py --outputDir STAR/SL1344

#Perform differential splicing analysis
cat splice_batches.txt | python ~/software/utils/submitJobs.py --MEM 5000 --jobname drimseq_differential_events --ncores 1 --queue normal --command "/software/R-3.3.0/bin/Rscript macrophage-gxe-study/SL1344/splicing/differentialEventExpression.R"


#Run FASTQTL on the event ratios
#Compress and index input files
bgzip results/SL1344/salmon/fastqtl_input/naive.norm_prop.txt && tabix -p bed results/SL1344/salmon/fastqtl_input/naive.norm_prop.txt.gz
bgzip results/SL1344/salmon/fastqtl_input/IFNg.norm_prop.txt && tabix -p bed results/SL1344/salmon/fastqtl_input/IFNg.norm_prop.txt.gz
bgzip results/SL1344/salmon/fastqtl_input/SL1344.norm_prop.txt && tabix -p bed results/SL1344/salmon/fastqtl_input/SL1344.norm_prop.txt.gz
bgzip results/SL1344/salmon/fastqtl_input/IFNg_SL1344.norm_prop.txt && tabix -p bed results/SL1344/salmon/fastqtl_input/IFNg_SL1344.norm_prop.txt.gz

#Run FastQTL on all conditions
cat results/SL1344/salmon/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname salmon_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --bed results/SL1344/salmon/fastqtl_input/naive.norm_prop.txt.gz --cov results/SL1344/salmon/fastqtl_input/naive.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon/fastqtl_output/naive_100kb_perm --execute True"
cat results/SL1344/salmon/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname salmon_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --bed results/SL1344/salmon/fastqtl_input/IFNg.norm_prop.txt.gz --cov results/SL1344/salmon/fastqtl_input/IFNg.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon/fastqtl_output/IFNg_100kb_perm --execute True"
cat results/SL1344/salmon/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --bed results/SL1344/salmon/fastqtl_input/SL1344.norm_prop.txt.gz --cov results/SL1344/salmon/fastqtl_input/SL1344.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon/fastqtl_output/SL1344_100kb_perm --execute True"
cat results/SL1344/salmon/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --bed results/SL1344/salmon/fastqtl_input/IFNg_SL1344.norm_prop.txt.gz --cov results/SL1344/salmon/fastqtl_input/IFNg_SL1344.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_perm --execute True"

#Merge chunks into single files
zcat results/SL1344/salmon/fastqtl_output/naive_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/naive_100kb_permuted.txt.gz
zcat results/SL1344/salmon/fastqtl_output/IFNg_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/IFNg_100kb_permuted.txt.gz
zcat results/SL1344/salmon/fastqtl_output/SL1344_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/SL1344_100kb_permuted.txt.gz
zcat results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_permuted.txt.gz

#Remove chunks
rm results/SL1344/salmon/fastqtl_output/*.chunk_*

#Get full p-values from fastQTL
cat results/SL1344/salmon/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --bed results/SL1344/salmon/fastqtl_input/naive.norm_prop.txt.gz --cov results/SL1344/salmon/fastqtl_input/naive.covariates_prop.txt --W 100000 --out results/SL1344/salmon/fastqtl_output/naive_100kb_full --execute True"
cat results/SL1344/salmon/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --bed results/SL1344/salmon/fastqtl_input/IFNg.norm_prop.txt.gz --cov results/SL1344/salmon/fastqtl_input/IFNg.covariates_prop.txt --W 100000 --out results/SL1344/salmon/fastqtl_output/IFNg_100kb_full --execute True"
cat results/SL1344/salmon/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --bed results/SL1344/salmon/fastqtl_input/SL1344.norm_prop.txt.gz --cov results/SL1344/salmon/fastqtl_input/SL1344.covariates_prop.txt --W 100000 --out results/SL1344/salmon/fastqtl_output/SL1344_100kb_full --execute True"
cat results/SL1344/salmon/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --bed results/SL1344/salmon/fastqtl_input/IFNg_SL1344.norm_prop.txt.gz --cov results/SL1344/salmon/fastqtl_input/IFNg_SL1344.covariates_prop.txt --W 100000 --out results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_full --execute True"

#Merge chunks into single files
zcat results/SL1344/salmon/fastqtl_output/naive_100kb_full.chunk_*.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/naive_100kb_pvalues.txt.gz
zcat results/SL1344/salmon/fastqtl_output/IFNg_100kb_full.chunk_*.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/IFNg_100kb_pvalues.txt.gz
zcat results/SL1344/salmon/fastqtl_output/SL1344_100kb_full.chunk_*.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/SL1344_100kb_pvalues.txt.gz
zcat results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_full.chunk_*.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_pvalues.txt.gz

#Remove chunks
rm results/SL1344/salmon/fastqtl_output/*.chunk_*

#Add SNP coordinates
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --fastqtl results/SL1344/salmon/fastqtl_output/naive_100kb_pvalues.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/naive_100kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --fastqtl results/SL1344/salmon/fastqtl_output/IFNg_100kb_pvalues.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/IFNg_100kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --fastqtl results/SL1344/salmon/fastqtl_output/SL1344_100kb_pvalues.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/SL1344_100kb_pvalues.coords.txt.gz"
echo hello | python ~/software/utils/submitJobs.py --MEM 5000 --jobname fastQTL_add_coords --command "python ~/software/utils/fastqtl/fastqtlAddSnpCoordinates.py --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --fastqtl results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_pvalues.txt.gz | bgzip > results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_pvalues.coords.txt.gz"

#Sort files by SNP coordinates
#awk command is necessary to change field separator from space to tab
zcat results/SL1344/salmon/fastqtl_output/naive_100kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/SL1344/salmon/fastqtl_output/naive_100kb_pvalues.sorted.txt.gz
zcat results/SL1344/salmon/fastqtl_output/IFNg_100kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/SL1344/salmon/fastqtl_output/IFNg_100kb_pvalues.sorted.txt.gz
zcat results/SL1344/salmon/fastqtl_output/SL1344_100kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/SL1344/salmon/fastqtl_output/SL1344_100kb_pvalues.sorted.txt.gz
zcat results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_pvalues.coords.txt.gz | awk -v OFS='\t' '{$1=$1; print $0}' | sort -k2,2 -k3,3n | bgzip > results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_pvalues.sorted.txt.gz
bsub -G team170 -n1 -R "span[hosts=1] select[mem>3000] rusage[mem=3000]" -q normal -M 3000 -n 3 -o FarmOut/sort_fastqtl.%J.jobout  "bash sort.sh"

#Index the output files using Tabix
tabix -s2 -b3 -e3 -f results/SL1344/salmon/fastqtl_output/naive_100kb_pvalues.coords.txt.gz
tabix -s2 -b3 -e3 -f results/SL1344/salmon/fastqtl_output/IFNg_100kb_pvalues.coords.txt.gz
tabix -s2 -b3 -e3 -f results/SL1344/salmon/fastqtl_output/SL1344_100kb_pvalues.coords.txt.gz
tabix -s2 -b3 -e3 -f results/SL1344/salmon/fastqtl_output/IFNg_SL1344_100kb_pvalues.coords.txt.gz


### Run FastQTL on Ensenbl transcript annotations
#Compress and index input files
bgzip results/SL1344/salmon_ensembl85/fastqtl_input/naive.norm_prop.txt && tabix -p bed results/SL1344/salmon_ensembl85/fastqtl_input/naive.norm_prop.txt.gz
bgzip results/SL1344/salmon_ensembl85/fastqtl_input/IFNg.norm_prop.txt && tabix -p bed results/SL1344/salmon_ensembl85/fastqtl_input/IFNg.norm_prop.txt.gz
bgzip results/SL1344/salmon_ensembl85/fastqtl_input/SL1344.norm_prop.txt && tabix -p bed results/SL1344/salmon_ensembl85/fastqtl_input/SL1344.norm_prop.txt.gz
bgzip results/SL1344/salmon_ensembl85/fastqtl_input/IFNg_SL1344.norm_prop.txt && tabix -p bed results/SL1344/salmon_ensembl85/fastqtl_input/IFNg_SL1344.norm_prop.txt.gz

#Run FastQTL on all conditions
cat results/SL1344/salmon_ensembl85/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname salmon_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --bed results/SL1344/salmon_ensembl85/fastqtl_input/naive.norm_prop.txt.gz --cov results/SL1344/salmon_ensembl85/fastqtl_input/naive.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon_ensembl85/fastqtl_output/naive_100kb_perm --execute True"
cat results/SL1344/salmon_ensembl85/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname salmon_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --bed results/SL1344/salmon_ensembl85/fastqtl_input/IFNg.norm_prop.txt.gz --cov results/SL1344/salmon_ensembl85/fastqtl_input/IFNg.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon_ensembl85/fastqtl_output/IFNg_100kb_perm --execute True"
cat results/SL1344/salmon_ensembl85/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --bed results/SL1344/salmon_ensembl85/fastqtl_input/SL1344.norm_prop.txt.gz --cov results/SL1344/salmon_ensembl85/fastqtl_input/SL1344.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon_ensembl85/fastqtl_output/SL1344_100kb_perm --execute True"
cat results/SL1344/salmon_ensembl85/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --bed results/SL1344/salmon_ensembl85/fastqtl_input/IFNg_SL1344.norm_prop.txt.gz --cov results/SL1344/salmon_ensembl85/fastqtl_input/IFNg_SL1344.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon_ensembl85/fastqtl_output/IFNg_SL1344_100kb_perm --execute True"

#Merge chunks into single files
zcat results/SL1344/salmon_ensembl85/fastqtl_output/naive_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon_ensembl85/fastqtl_output/naive_100kb_permuted.txt.gz
zcat results/SL1344/salmon_ensembl85/fastqtl_output/IFNg_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon_ensembl85/fastqtl_output/IFNg_100kb_permuted.txt.gz
zcat results/SL1344/salmon_ensembl85/fastqtl_output/SL1344_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon_ensembl85/fastqtl_output/SL1344_100kb_permuted.txt.gz
zcat results/SL1344/salmon_ensembl85/fastqtl_output/IFNg_SL1344_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon_ensembl85/fastqtl_output/IFNg_SL1344_100kb_permuted.txt.gz

#Remove chunks
rm results/SL1344/salmon_ensembl85/fastqtl_output/*.chunk_*



#Run FastQTL on full Ensembl transcript annotations without filtering
#Compress and index input files
bgzip results/SL1344/salmon_ensembl85_full/fastqtl_input/naive.norm_prop.txt && tabix -p bed results/SL1344/salmon_ensembl85_full/fastqtl_input/naive.norm_prop.txt.gz
bgzip results/SL1344/salmon_ensembl85_full/fastqtl_input/IFNg.norm_prop.txt && tabix -p bed results/SL1344/salmon_ensembl85_full/fastqtl_input/IFNg.norm_prop.txt.gz
bgzip results/SL1344/salmon_ensembl85_full/fastqtl_input/SL1344.norm_prop.txt && tabix -p bed results/SL1344/salmon_ensembl85_full/fastqtl_input/SL1344.norm_prop.txt.gz
bgzip results/SL1344/salmon_ensembl85_full/fastqtl_input/IFNg_SL1344.norm_prop.txt && tabix -p bed results/SL1344/salmon_ensembl85_full/fastqtl_input/IFNg_SL1344.norm_prop.txt.gz

#Run FastQTL on all conditions
cat results/SL1344/salmon_ensembl85_full/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname salmon_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --bed results/SL1344/salmon_ensembl85_full/fastqtl_input/naive.norm_prop.txt.gz --cov results/SL1344/salmon_ensembl85_full/fastqtl_input/naive.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon_ensembl85_full/fastqtl_output/naive_100kb_perm --execute True"
cat results/SL1344/salmon_ensembl85_full/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname salmon_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg.ASE.vcf.gz --bed results/SL1344/salmon_ensembl85_full/fastqtl_input/IFNg.norm_prop.txt.gz --cov results/SL1344/salmon_ensembl85_full/fastqtl_input/IFNg.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon_ensembl85_full/fastqtl_output/IFNg_100kb_perm --execute True"
cat results/SL1344/salmon_ensembl85_full/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/SL1344.ASE.vcf.gz --bed results/SL1344/salmon_ensembl85_full/fastqtl_input/SL1344.norm_prop.txt.gz --cov results/SL1344/salmon_ensembl85_full/fastqtl_input/SL1344.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon_ensembl85_full/fastqtl_output/SL1344_100kb_perm --execute True"
cat results/SL1344/salmon_ensembl85_full/fastqtl_input/all_chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname leafcutter_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/IFNg_SL1344.ASE.vcf.gz --bed results/SL1344/salmon_ensembl85_full/fastqtl_input/IFNg_SL1344.norm_prop.txt.gz --cov results/SL1344/salmon_ensembl85_full/fastqtl_input/IFNg_SL1344.covariates_prop.txt --W 100000 --permute '100 10000' --out results/SL1344/salmon_ensembl85_full/fastqtl_output/IFNg_SL1344_100kb_perm --execute True"

#Merge chunks into single files
zcat results/SL1344/salmon_ensembl85_full/fastqtl_output/naive_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon_ensembl85_full/fastqtl_output/naive_100kb_permuted.txt.gz
zcat results/SL1344/salmon_ensembl85_full/fastqtl_output/IFNg_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon_ensembl85_full/fastqtl_output/IFNg_100kb_permuted.txt.gz
zcat results/SL1344/salmon_ensembl85_full/fastqtl_output/SL1344_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon_ensembl85_full/fastqtl_output/SL1344_100kb_permuted.txt.gz
zcat results/SL1344/salmon_ensembl85_full/fastqtl_output/IFNg_SL1344_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/salmon_ensembl85_full/fastqtl_output/IFNg_SL1344_100kb_permuted.txt.gz

#Remove chunks
rm results/SL1344/salmon_ensembl85_full/fastqtl_output/*.chunk_*

