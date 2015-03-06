#Process genotyping data
cut -f2 rna_seq_genptypes.txt | tail -n +2 > rna_seq_genotype_names.txt

#Only keep the relevant genotypes from the vcf file
bcftools view -O z -S rna_seq_genotype_names.txt hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.vcf.gz > hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.subset.vcf.gz 

#Only keep SNPs with MAF > 0.05
bcftools filter -O z -i 'MAF[0] >= 0.05' hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.subset.vcf.gz > hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.subset.maf005.vcf.gz 

#Convert vcf to bcf
bcftools view -O b hipsci.wec.gtarray.HumanCoreExome-12_v1_0.858samples.20141111.genotypes.subset.maf005.vcf.gz > selected_genotypes.gz

#Change the chromosome names of the bcf file
bcftools annotate --rename-chrs GRCh38ToHg38_chromosome_map.txt -Ob selected_genotypes.bcf.gz > selected_genotypes.hg19.bcf.gz 

#Construct dictionary of the hg38 ref sequence
java -Xmx1000m -jar ~/software/picard-tools-1.129/picard.jar CreateSequenceDictionary R= hg38.fa O= hg38.dict
java -Xmx1000m -jar ~/software/picard-tools-1.129/picard.jar CreateSequenceDictionary R= hg19.fa O= hg19.dict

#Convert BCF to VCF
bcftools view selected_genotypes.hg19.bcf.gz > selected_genotypes.hg19.vcf

#Lift over variants from hg19 to hg38
/software/java/bin/java -Xmx1000m -jar ~/software/GenomeAnalysisTK.jar -T LiftoverVariants --newSequenceDictionary ../../../annotations/hg38/hg38.fa --chain hg19ToHg38.over.chain --out selected_genotypes.hg38.vcf --variant selected_genotypes.hg19.reordered.vcf -R ../../../annotations/hg19/hg19.fa 