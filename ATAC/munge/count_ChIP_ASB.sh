#Analyse PU.1 data from Schulze
cut -f1 macrophage-gxe-study/data/chromatin/ASB/Schultze_2016.txt | head -n3 | python ~/software/utils/submitJobs.py --MEM 6000 --jobname bamCountASE --command "python ~/software/rasqual/scripts/bamCountASE.py --indir processed/Schultze_2016 --outdir processed/Schultze_2016 --insuffix .no_duplicates.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True --Xmx 4g"

#Analyse CEBPb data from O'Callaghan
cut -f1 macrophage-gxe-study/data/chromatin/ASB/OCallaghan.txt | tail -n 7| python ~/software/utils/submitJobs.py --MEM 6000 --jobname bamCountASE --command "python ~/software/rasqual/scripts/bamCountASE.py --indir processed/OCallaghan --outdir processed/OCallaghan --insuffix .no_duplicates.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True --Xmx 4g"

#Aanlyse PU.1 and CEBPb from Rehli
cut -f1 macrophage-gxe-study/data/chromatin/ASB/Rehli.txt | python ~/software/utils/submitJobs.py --MEM 6000 --jobname bamCountASE --command "python ~/software/rasqual/scripts/bamCountASE.py --indir processed/Rehli --outdir processed/Rehli --insuffix .no_duplicates.bam --reference ../../annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sites genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz --execute True --Xmx 4g"