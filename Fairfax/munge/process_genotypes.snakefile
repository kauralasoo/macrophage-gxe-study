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
		protected("processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz")
	resources:
		mem = 1000
	threads: 3
	shell:
		"~/software/vcflib/bin/vcfuniq {input} | bcftools norm -m+any - | "
		"bcftools view -m2 -M2 - | "
		"bcftools annotate -O z --set-id +'%CHROM\_%POS' > {output}"

rule extract_variant_infromation:
	input:
		"processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz"
	output:
		"processed/Fairfax/merged_genotypes/fairfax_genotypes.variant_information.txt.gz"
	resources:
		mem = 1000
	threads: 2	
	shell:
		"bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\n' {input} | bgzip > {output}"

rule index_full_vcf:
	input:
		"processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz"
	output:
		"processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz.csi"
	resources:
		mem = 1000
	threads: 1
	shell:
		"bcftools index {input}"

rule extract_chr:
	input:
		vcf = "processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz",
		index = "processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz.csi"
	output:
		vcf = "processed/Fairfax/geno_by_chr/{chr}.vcf.gz"
	resources:
		mem = 1000
	threads: 1
	shell:
		"bcftools view -O z -r {wildcards.chr} {input.vcf} > {output.vcf}"

rule convert_vcf_to_gds:
	input:
		vcf = "processed/Fairfax/geno_by_chr/{chr}.vcf.gz"
	output:
		gds = "processed/Fairfax/geno_by_chr/{chr}.gds"
	resources:
		mem = 1000
	threads: 1
	shell:
		"/software/R-3.4.0/bin/Rscript ~/software/utils/vcf/vcfToGds.R -v {input.vcf} -g {output.gds}"


rule make_all:
	input:
		"processed/Fairfax/merged_genotypes/fairfax_genotypes.variant_information.txt.gz",
		expand("processed/Fairfax/geno_by_chr/{chr}.vcf.gz", chr=chromosomes),
		expand("processed/Fairfax/geno_by_chr/{chr}.gds", chr=chromosomes)
	output:
		"processed/Fairfax/out.txt"
	resources:
		mem = 1000
	threads: 1
	shell:
		"echo 'Done!' > {output}"



