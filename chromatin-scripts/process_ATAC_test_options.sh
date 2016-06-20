#test-run rasqual
echo "ATAC_peak_55687" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/IFNg.expression.bin --offsets results/ATAC/rasqual/input/IFNg.offsets.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.no_quote.txt --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_counts.txt --execute True > CTSC_enhancher.cond_A.rasqual"
python ~/software/utils/rasqualToIGV.py --rasqualOut CTSC_enhancher.cond_A.rasqual --geneid ATAC_peak_55687 > CTSC_enhancer.cond_A.gwas


echo "ATAC_peak_55687" | python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.offsets.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.no_quote.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_counts.txt --execute True > CTSC_enhancher.cond_A.rasqual

echo -e "batch_1\tATAC_peak_55686,ATAC_peak_55687" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/IFNg.expression.bin --offsets results/ATAC/rasqual/input/IFNg.offsets.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_counts.txt --outdir results/ATAC/rasqual/output/ --execute True"

echo -e "batch_1_2kb\tATAC_peak_55686,ATAC_peak_55687" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/IFNg.expression.bin --offsets results/ATAC/rasqual/input/IFNg.offsets.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/IFNg.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_2kb.txt --outdir results/ATAC/rasqual/output/ --execute True"

cat results/ATAC/rasqual/input/peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.offsets.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_50kb.txt --outdir results/ATAC/rasqual/output_50kb/ --execute True"


cat results/ATAC/rasqual/output_50kb/batch_*.txt | cut -f 1,2,3,4,11,12,17,18,23 | bgzip > results/ATAC/rasqual/output_50kb/naive_50kb.txt.gz

#Use library size instead of RLE
cat results/ATAC/rasqual/input/peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.library_size.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_50kb.txt --outprefix results/ATAC/rasqual/output/naive_50kb_ls --execute True"

cat results/ATAC/rasqual/output/naive_50kb_ls.batch_*.txt | cut -f 1,2,3,4,11,12,17,18,23 | bgzip > results/ATAC/rasqual/output/naive_50kb_ls.txt.gz

#Use 5 PEER covariates
cat results/ATAC/rasqual/input/peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.library_size.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_50kb.txt --covariates results/ATAC/rasqual/input/naive.covariates5.bin --outprefix results/ATAC/rasqual/output/naive_50kb_cov --execute True"

cat results/ATAC/rasqual/output/naive_50kb_cov.batch_*.txt | cut -f 1,2,3,4,11,12,17,18,23 | bgzip > results/ATAC/rasqual/output/naive_50kb_ls.txt.gz


#Run rasqual on chr11
cat results/ATAC/rasqual/input/peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.library_size.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_50kb.txt --outprefix results/ATAC/rasqual/output/naive_chr11_50kb --execute True"
cat results/ATAC/rasqual/output/naive_chr11_50kb.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 | bgzip > results/ATAC/rasqual/output/naive_chr11_50kb.txt.gz

#Run rasqual on chr11 (500kb window)
cat results/ATAC/rasqual/input/peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.library_size.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_500kb.txt --outprefix results/ATAC/rasqual/output/naive_chr11_500kb --execute True"
cat results/ATAC/rasqual/output/naive_chr11_500kb.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 | bgzip > results/ATAC/rasqual/output/naive_chr11_500kb.txt.gz


#Run rasqual on chr11 (GC-corrected)
cat results/ATAC/rasqual/input/peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.lib_size_gc.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_50kb.txt --outprefix results/ATAC/rasqual/output/naive_chr11_50kb_gc --execute True"
cat results/ATAC/rasqual/output/naive_chr11_50kb_gc.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 | bgzip > results/ATAC/rasqual/output/naive_chr11_50kb_gc.txt.gz


#Run rasqual on chr11 (with 2 covariates)
cat results/ATAC/rasqual/input/peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.library_size.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_50kb.txt --outprefix results/ATAC/rasqual/output/naive_chr11_50kb_cov --covariates results/ATAC/rasqual/input/naive.covariates_pc2.bin --execute True"
cat results/ATAC/rasqual/output/naive_chr11_50kb_cov.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 | bgzip > results/ATAC/rasqual/output/naive_chr11_50kb_cov.txt.gz

#Run rasqual on chr11 (with 7 covariates)
cat results/ATAC/rasqual/input/peak_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual --command "python ~/software/utils/runRasqual.py --readCounts results/ATAC/rasqual/input/naive.expression.bin --offsets results/ATAC/rasqual/input/naive.library_size.bin --n 27 --geneids results/ATAC/rasqual/input/peak_names.txt --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/ATAC/rasqual/input/peak_snp_count_50kb.txt --outprefix results/ATAC/rasqual/output/naive_chr11_50kb_cov_pc4 --covariates results/ATAC/rasqual/input/naive.covariates_pc4.bin --execute True"
cat results/ATAC/rasqual/output/naive_chr11_50kb_cov_pc4.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 | bgzip > results/ATAC/rasqual/output/naive_chr11_50kb_cov_pc4.txt.gz


#Count fragment lengths
cut -f1 macrophage-chromatin/data/SL1344/ATAC_Salmonella_names.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname bedCountFragmentLengths --command "python ~/software/utils/coverage/bedCountFragmentLengths.py  --indir processed/SL1344 --outdir processed/SL1344 --execute True"

#Calculate GC content for ATAC peaks
bedtools nuc -fi ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed annotations/ATAC_Seq_joint_peaks.gff3 > annotations/ATAC_Seq_joint_peaks.nuc_content.txt
cut -f9,11 annotations/ATAC_Seq_joint_peaks.nuc_content.txt | cut -f2 -d= > annotations/ATAC_Seq_joint_peaks.GC_content.txt

#Rsync small files to a different directory for copying
rsync -R SL1344/*/*.fragment_lengths.txt SL1344_rsync/
rsync -R SL1344/*/*.MarkDuplicates.txt SL1344_rsync/
rsync -R SL1344/*/*.chr_counts SL1344_rsync/
rsync -R SL1344/*/*.verifyBamID.bestSM SL1344_rsync/
rsync -R SL1344/*/*.consensus_peaks.counts.* SL1344_rsync/



### FASTQTL ###

#Zip+index bed files
bgzip results/ATAC/fastqtl/input/naive.expression.txt && tabix -p bed results/ATAC/fastqtl/input/naive.expression.txt.gz
bgzip results/ATAC/fastqtl/input/IFNg.expression.txt && tabix -p bed results/ATAC/fastqtl/input/IFNg.expression.txt.gz
bgzip results/ATAC/fastqtl/input/SL1344.expression.txt && tabix -p bed results/ATAC/fastqtl/input/SL1344.expression.txt.gz
bgzip results/ATAC/fastqtl/input/IFNg_SL1344.expression.txt && tabix -p bed results/ATAC/fastqtl/input/IFNg_SL1344.expression.txt.gz

#Run FastQTL with CQN expression values and without covariates
cat results/ATAC/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression.txt.gz --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output/naive_chr11_50kb_perm --execute True"
zcat results/ATAC/fastqtl/output/naive_chr11_50kb_perm.chunk_*.txt.gz |  bgzip > results/ATAC/fastqtl/output/naive_chr11_50kb_perm.txt.gz

cat results/ATAC/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 3000 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression.txt.gz --W 50000 --out results/ATAC/fastqtl/output/naive_50kb_full --execute True"
zcat results/ATAC/fastqtl/output/naive_50kb_full.chunk_*.txt.gz | bgzip > results/ATAC/fastqtl/output/naive_50kb_full.txt.gz

#Run FastQTL wtih TPM expression values
bgzip results/ATAC/fastqtl/input/naive.expression_tpm.txt && tabix -p bed results/ATAC/fastqtl/input/naive.expression_tpm.txt.gz
bgzip results/ATAC/fastqtl/input/IFNg.expression_tpm.txt && tabix -p bed results/ATAC/fastqtl/input/IFNg.expression_tpm.txt.gz
bgzip results/ATAC/fastqtl/input/SL1344.expression_tpm.txt && tabix -p bed results/ATAC/fastqtl/input/SL1344.expression_tpm.txt.gz
bgzip results/ATAC/fastqtl/input/IFNg_SL1344.expression_tpm.txt && tabix -p bed results/ATAC/fastqtl/input/IFNg_SL1344.expression_tpm.txt.gz

cat results/ATAC/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 200 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression_tpm.txt.gz --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output/naive_chr11_50kb_tpm_perm --execute True"
zcat results/ATAC/fastqtl/output/naive_chr11_50kb_tpm_perm.chunk_*.txt.gz |  bgzip > results/ATAC/fastqtl/output/naive_chr11_50kb_tpm_perm.txt.gz
rm results/ATAC/fastqtl/output/*chunk_*

#Run Fastqtl with 7 covariates covariates (4 PCs)
cat results/ATAC/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression_tpm.txt.gz --cov results/ATAC/fastqtl/input/naive.covariates_pc4.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output/naive_chr11_50kb_tpm_perm_pc4 --execute True"
zcat results/ATAC/fastqtl/output/naive_chr11_50kb_tpm_perm_pc4.chunk_*.txt.gz |  bgzip > results/ATAC/fastqtl/output/naive_chr11_50kb_tpm_perm_pc4.txt.gz
rm results/ATAC/fastqtl/output/*chunk_*

#Run fastQTL with 3 covariates (no PCs)
cat results/ATAC/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression_tpm.txt.gz --cov results/ATAC/fastqtl/input/naive.covariates_pc0.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output/naive_chr11_50kb_tpm_perm_pc0 --execute True"
zcat results/ATAC/fastqtl/output/naive_chr11_50kb_tpm_perm_pc0.chunk_*.txt.gz |  bgzip > results/ATAC/fastqtl/output/naive_chr11_50kb_tpm_perm_pc0.txt.gz
rm results/ATAC/fastqtl/output/*chunk_*

#Run Fastqtl with 3 peer covariates
cat results/ATAC/fastqtl/input/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/ATAC/rasqual/input/naive.ASE.vcf.gz --bed results/ATAC/fastqtl/input/naive.expression.txt.gz --cov results/ATAC/fastqtl/input/naive.covariates_peer3.txt --W 50000 --permute '100 10000' --out results/ATAC/fastqtl/output/naive_peer3_50kb_perm --execute True"
zcat results/ATAC/fastqtl/output/naive_peer3_50kb_perm.chunk_*.txt.gz |  bgzip > results/ATAC/fastqtl/output/naive_peer3_50kb_permuted.txt.gz
rm results/ATAC/fastqtl/output/*chunk_*