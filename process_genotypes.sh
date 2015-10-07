#Process genotyping data
cut -f2 rna_seq_genptypes.txt | tail -n +2 > rna_seq_genotype_names.txt

#Only keep the relevant samples from the vcf file
bcftools view -O z -S rna_seq_genotype_names.txt hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.vcf.gz > hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.subset.vcf.gz 

#Only keep SNPs with MAF > 0.05
bcftools filter -i 'MAF[0] >= 0.05' hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.subset.vcf.gz > selected_genotypes.vcf

#Change the chromosome names of the bcf file
bcftools annotate --rename-chrs GRCh38ToHg38_chromosome_map.txt selected_genotypes.vcf > selected_genotypes.hg19.vcf 

#Lift over variants from hg19 to hg38 (Only works on my MacBook Air)
CrossMap.py vcf hg19ToHg38.over.chain selected_genotypes.hg19.vcf hg38.fa selected_genotypes.hg38.vcf

#Fix chromosome names and remove ALT contigs introduced by CrossMap
python ../macrophage-gxe-study/genotypes/postprocessCrossmap.py --vcf selected_genotypes.hg38.vcf > selected_genotypes.hg38.fixed.vcf

#Convert chromosome names back to GRCh38
bcftools annotate --rename-chrs Hg38ToGRCh38_chromosome_map.txt selected_genotypes.hg38.fixed.vcf > selected_genotypes.GRCh38.vcf 

#Sort vcf file by position
grep '^#' selected_genotypes.GRCh38.vcf > selected_genotypes.GRCh38.sorted.vcf && grep -v '^#' selected_genotypes.GRCh38.vcf | LC_ALL=C sort -k1,1 -k2,2n >> selected_genotypes.GRCh38.sorted.vcf

#Compress and index the VCF file
bcftools view -O z selected_genotypes.GRCh38.sorted.vcf > selected_genotypes.GRCh38.sorted.vcf.gz
bcftools index selected_genotypes.GRCh38.sorted.vcf.gz

#Find the genotyped SNPs from the imputed VCF file (will get rid of NAs)
bcftools view -R array_genotypes.59_samples.vcf.gz -o array_genotypes.59_samples.imputed.vcf.gz -O z autosomes.snps_only.no_replicates.vcf.gz 

#Keep only unique genotypes
python ~/software/utils/vcf/vcfDetectDuplicateVariants.py --vcf array_genotypes.59_samples.imputed.vcf.gz --duplicates dup.txt | bgzip > array_genotypes.59_samples.imputed.uniq.vcf.gz 