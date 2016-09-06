#Construct transcriptome index for Salmon
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o construct_index.%J.jobout "salmon --no-version-check index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i salmon_index_85"

#Quantify Ensembl transcript expression using Salmon
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/submitJobs.py --MEM 10000 --jobname salmon_quant --ncores 8 --queue normal --command "python ~/software/utils/align/salmonQuant.py --outputDir processed/SL1344_salmon/ --outputSubdir ensembl_full_bootstraps_50 --fastqDir fastq/SL1344 --index ../../annotations/GRCh38/genes/Ensembl_85/salmon_index_85/ --libType ISR --geneMap ../../annotations/GRCh38/genes/Ensembl_85/Homo_sapiens.GRCh38.85.gtf --nCores 8 --numBootstraps 50"

#Copy gene and transcript estimates back to the main results folder
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/copySalmonOutput.py --currentDir processed/SL1344_salmon/ --currentSubdir ensembl_full_bootstraps_50 --newDir STAR/SL1344/ --filename quant.sf --suffix ensembl85
cut -f1 macrophage-gxe-study/data/sample_lists/SL1344/SL1344_names_all.txt | python ~/software/utils/copySalmonOutput.py --currentDir processed/SL1344_salmon/ --currentSubdir ensembl_full_bootstraps_50 --newDir STAR/SL1344/ --filename quant.genes.sf --suffix ensembl85