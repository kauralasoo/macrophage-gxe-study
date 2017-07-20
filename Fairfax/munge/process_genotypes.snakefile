chromosomes = [10,11,12,13,14,15,16,17,18,19,1,20,21,22,2,3,4,5,6,7,8,9]

rule filter_maf_info:
	input:
		vcf = "processed/Fairfax/raw_genotypes/{chr}.pbwt_reference_impute.vcf.gz"
	output:
		vcf = temp("processed/Fairfax/filtered_genotypes/{chr}.maf_info.vcf.gz")
	resources:
		mem = 1000
	threads: 1
	shell:
		"bcftools filter -i 'INFO[0] >= 0.4 & MAF[0] >= 0.05' {input.vcf} > {output.vcf}"

rule filter_hwe:
	input:
		vcf = "processed/Fairfax/filtered_genotypes/{chr}.maf_info.vcf.gz"
	output:
		vcf = temp("processed/Fairfax/filtered_genotypes/{chr}.maf_info_hwe.vcf.gz")
	resources:
		mem = 1000
	threads: 2
	shell:
		"vcftools --gzvcf {input.vcf} --hwe 1e-6 --recode --recode-INFO-all --stdout | bcftools view -O z - > {output.vcf}"

rule index_vcf:
	input:
		vcf = "processed/Fairfax/filtered_genotypes/{chr}.maf_info_hwe.vcf.gz"
	output:
		temp("processed/Fairfax/filtered_genotypes/{chr}.maf_info_hwe.vcf.gz.csi")
	resources:
		mem = 1000
	threads: 1
	shell:
		"bcftools index {input.vcf}"


rule concat_chromosomes:
	input:
		vcfs = expand("processed/Fairfax/filtered_genotypes/{chr}.maf_info_hwe.vcf.gz", chr=chromosomes),
		indices = expand("processed/Fairfax/filtered_genotypes/{chr}.maf_info_hwe.vcf.gz.csi", chr=chromosomes)
	output:
		temp("processed/Fairfax/merged_genotypes/fairfax_genotypes.vcf")
	resources:
		mem = 1000
	threads: 1
	shell:
		"bcftools concat -a {input.vcfs} > {output}"

rule sort_vcf:
	input:
		"processed/Fairfax/merged_genotypes/fairfax_genotypes.vcf"
	output:
		temp("processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.vcf")
	resources:
		mem = 1000
	threads: 1
	shell:
		"~/software/vcflib/bin/vcfsort {input} > {output}"

rule keep_unique_biallelic:
	input:
		"processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.vcf"
	output:
		"processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz"
	resources:
		mem = 1000
	threads: 3
	shell:
		"~/software/vcflib/bin/vcfuniq {input} | bcftools norm -m+any - | "
		"bcftools view -m2 -M2 - | "
		"bcftools annotate -O z --set-id +'%CHROM\_%POS' > {output}"





