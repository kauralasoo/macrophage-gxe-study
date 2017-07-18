chromosomes = [10,11,12,13,14,15,16,17,18,19,1,20,21,22,2,3,4,5,6,7,8,9]

filter_maf_info:
	input:
		vcf = "processed/raw_genotypes/{chr}.pbwt_reference_impute.vcf.gz"
	output:
		vcf = "processed/filtered_genotypes/{chr}.maf_info.vcf.gz"
	resources:
		mem = 1000
	threads: 1
	shell:
		"bcftools filter -i 'INFO[0] >= 0.4 & MAF[0] >= 0.05' {input.vcf} |" 
		"bcftools annotate -O z --set-id +'%CHROM\_%POS' > {output.vcf}"

filter_hwe:
	input:
		vcf = "processed/filtered_genotypes/{chr}.maf_info.vcf.gz"
	output:
		vcf = "processed/filtered_genotypes/{chr}.maf_info_hwe.vcf.gz"
	resources:
		mem = 1000
	threads: 1
	shell:
		"vcftools --gzvcf {input.vcf} --hwe 1e-6 --recode --recode-INFO-all --stdout | bcftools view -O z - > {output.vcf}"

index_vcf:
	input:
		vcf = "processed/filtered_genotypes/{chr}.maf_info_hwe.vcf.gz"
	output:
		"processed/filtered_genotypes/{chr}.maf_info_hwe.vcf.gz.csi"
	resources:
		mem = 1000
	threads: 1
	shell:
		"bcftools index {input.vcf}"


concat_chromosomes:
	input:
		expand("processed/filtered_genotypes/{chr}.maf_info_hwe.vcf.gz", chr=chromosomes)
	output:
		"processed/merged/fairfax_genotypes.vcf"
	resources:
		mem = 1000
	threads: 1
	shell:
		"bcftools concat -a {input} > {output}"




