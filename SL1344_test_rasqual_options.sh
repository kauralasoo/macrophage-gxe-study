bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o eigenMT.%J.jobout "python eigenMT.py --CHROM 19 --QTL qtls.txt --GEN genotypes.txt --GENPOS gen.positions.txt --PHEPOS phe.positions.txt --OUT exampleOut.txt"

#Plain RASQUAL (no covariates, library sizes), 100kb
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb --execute True"

bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 -o FarmOut/rasqual_concat_results.%J.jobout "cat results/SL1344/rasqual/output/chr11_naive_100kb.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 > results/SL1344/rasqual/output/chr11_naive_100kb.txt"

#RASQUAL (no covariates, library sizes + GC), 100kb
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb_gc --execute True"

bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 -o FarmOut/rasqual_concat_results.%J.jobout "cat results/SL1344/rasqual/output/chr11_naive_100kb_gc.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 > results/SL1344/rasqual/output/chr11_naive_100kb_gc.txt"

#RASQUAL (SVD covariates, library sizes + GC), 100kb
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_svd --covariates results/SL1344/rasqual/input/naive.svd_covariates.bin --execute True"

bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 -o FarmOut/rasqual_concat_results.%J.jobout "cat results/SL1344/rasqual/output/chr11_naive_100kb_gc_svd.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 > results/SL1344/rasqual/output/chr11_naive_100kb_gc_svd.txt"

#Plain RASQUAL (no covariates, RLE sizes), 100kb
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.RLE_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb_RLE --execute True"

bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 -o FarmOut/rasqual_concat_results.%J.jobout "cat results/SL1344/rasqual/output/chr11_naive_100kb_RLE.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 > results/SL1344/rasqual/output/chr11_naive_100kb_RLE.txt"

#Plain RASQUAL (no covariates, RLE sizes + GC), 100kb
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.gc_RLE_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_RLE --execute True"

bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 -o FarmOut/rasqual_concat_results.%J.jobout "cat results/SL1344/rasqual/output/chr11_naive_100kb_gc_RLE.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 > results/SL1344/rasqual/output/chr11_naive_100kb_gc_RLE.txt"

#RASQUAL (4 PEER covariates + sex, library sizes + GC), 100kb
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER --covariates results/SL1344/rasqual/input/naive.PEER_covariates.bin --execute True"

bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 -o FarmOut/rasqual_concat_results.%J.jobout "cat results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 > results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER.txt"

#RASQUAL (4 PEER covariates + sex, library sizes + GC), 100kb (Using new VCF file with updated snp_ids and new rasqual binary)
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_renamed --covariates results/SL1344/rasqual/input/naive.PEER_covariates.bin --rasqualBin rasqual --execute True"

echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/utils/rasqual/mergeRasqualBatches.py --prefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_renamed"

#RASQUAL (2 PEER covariates + sex, library sizes + GC), 100kb (Using new VCF file with updated snp_ids and new rasqual binary)
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_n3 --covariates results/SL1344/rasqual/input/naive.PEER_covariates_n3.bin --rasqualBin rasqual --execute True"

echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/utils/rasqual/mergeRasqualBatches.py --prefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_n3"


#RASQUAL (population only, 2 PEER covariates + sex, library sizes + GC), 100kb (Using new VCF file with updated snp_ids and new rasqual binary)
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_n3_pop --covariates results/SL1344/rasqual/input/naive.PEER_covariates_n3.bin --rasqualBin rasqual --parameters '\--population-only' --execute True"

echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/utils/rasqual/mergeRasqualBatches.py --prefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_n3_pop"

#RASQUAL (population only, library sizes + GC), 100kb (Using new VCF file with updated snp_ids and new rasqual binary)
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_pop --rasqualBin rasqual --parameters '\--population-only' --execute True"

echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/utils/rasqual/mergeRasqualBatches.py --prefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_pop"



#RASQUAL (population only, SVD covariates, library sizes + GC), 100kb (Using new VCF file with updated snp_ids and new rasqual binary)
cat results/SL1344/rasqual/input/chr11_batches.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname runRasqual_all --command "python ~/software/utils/rasqual/runRasqual.py --readCounts results/SL1344/rasqual/input/naive.expression.bin  --offsets results/SL1344/rasqual/input/naive.gc_library_size.bin --n 69 --geneids results/SL1344/rasqual/input/feature_names.txt --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --geneMetadata results/SL1344/rasqual/input/gene_snp_count_100kb.txt --outprefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_svd_pop --covariates results/SL1344/rasqual/input/naive.svd_covariates.bin --rasqualBin rasqual --parameters '\--population-only' --execute True"

echo "merge" | python ~/software/utils/submitJobs.py --MEM 1000 --jobname mergeRasqualBatches --command "python ~/software/utils/rasqual/mergeRasqualBatches.py --prefix results/SL1344/rasqual/output/chr11_naive_100kb_gc_svd_pop"




bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o eigenMT.%J.jobout "python ~/software/eigenMTwithTestData/eigenMT.py --CHROM 11 --QTL qtls_2.txt --GEN genotypes.txt --GENPOS snp_positions.txt --PHEPOS gene_positions.txt --OUT exampleOut.txt --cis_dist 1e5 --external"



#Compare to FastQTL
bgzip results/SL1344/fastqtl/input_chr11/naive.cqn_expression.txt && tabix -p bed results/SL1344/fastqtl/input_chr11/naive.cqn_expression.txt.gz
bgzip results/SL1344/fastqtl/input_chr11/IFNg.cqn_expression.txt && tabix -p bed results/SL1344/fastqtl/input_chr11/IFNg.cqn_expression.txt.gz
bgzip results/SL1344/fastqtl/input_chr11/SL1344.cqn_expression.txt && tabix -p bed results/SL1344/fastqtl/input_chr11/SL1344.cqn_expression.txt.gz
bgzip results/SL1344/fastqtl/input_chr11/IFNg_SL1344.cqn_expression.txt && tabix -p bed results/SL1344/fastqtl/input_chr11/IFNg_SL1344.cqn_expression.txt.gz

#Run FastQTL without covariates
cat results/SL1344/fastqtl/input_chr11/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --bed results/SL1344/fastqtl/input_chr11/naive.cqn_expression.txt.gz --W 100000 --permute '100 10000' --out results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_perm --execute True"
zcat results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_perm.txt.gz

cat results/SL1344/fastqtl/input_chr11/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --bed results/SL1344/fastqtl/input_chr11/naive.cqn_expression.txt.gz --W 100000 --out results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_full --execute True"
zcat results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_full.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_full.txt.gz

#Run FastQTL with covariates (sex + 4xPEER)
cat results/SL1344/fastqtl/input_chr11/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --bed results/SL1344/fastqtl/input_chr11/naive.cqn_expression.txt.gz --W 100000 --permute '100 10000' --cov results/SL1344/fastqtl/input_chr11/naive.covariates.txt  --out results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_cov_perm --execute True"
zcat results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_cov_perm.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_cov_perm.txt.gz

cat results/SL1344/fastqtl/input_chr11/chunk_table.txt | python ~/software/utils/submitJobs.py --MEM 500 --jobname run_fastQTL --ncores 1 --command "python ~/software/utils/fastqtl/runFastQTL.py --vcf results/SL1344/rasqual/input/naive.ASE.vcf.gz --bed results/SL1344/fastqtl/input_chr11/naive.cqn_expression.txt.gz --W 100000 --cov results/SL1344/fastqtl/input_chr11/naive.covariates.txt  --out results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_cov_full --execute True"
zcat results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_cov_full.chunk_*.txt.gz | bgzip > results/SL1344/fastqtl/output_chr11/naive_chr11_100kb_cov_full.txt.gz





