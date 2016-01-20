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

bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -M 1000 -o FarmOut/rasqual_concat_results.%J.jobout "cat results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_renamed.batch_*.txt | cut -f 1,2,3,4,7,11,12,17,18,23 > results/SL1344/rasqual/output/chr11_naive_100kb_gc_PEER_renamed.txt"




bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o eigenMT.%J.jobout "python ~/software/eigenMTwithTestData/eigenMT.py --CHROM 11 --QTL qtls_2.txt --GEN genotypes.txt --GENPOS snp_positions.txt --PHEPOS gene_positions.txt --OUT exampleOut.txt --cis_dist 1e5 --external"


