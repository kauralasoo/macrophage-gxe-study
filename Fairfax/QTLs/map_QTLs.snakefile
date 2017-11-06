rule map_qtls:
	input:
		expand("processed/{{study}}/qtltools/output/{annot_type}/{condition}.permuted.txt.gz", annot_type = config["annot_type"], condition = config["conditions"]),
		expand("processed/{{study}}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz", annot_type = config["annot_type"], condition = config["conditions"]),
		expand("processed/{{study}}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz.tbi", annot_type = config["annot_type"], condition = config["conditions"]),
		expand("processed/{{study}}/eigenMT/output/{condition}.{chromosome}.txt", condition = config["conditions"], chromosome = config["chromosomes"])
	output:
		"processed/{study}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"

#Compres and index input bed file
rule compress_bed:
	input:
		bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt"
	output:
		bed = protected("processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz"),
		bed_index = protected("processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz.tbi")
	threads: 1
	resources:
		mem = 100
	shell:
		"bgzip {input.bed} && tabix -p bed {output.bed}"

#Run QTLtools in permutation mode
rule permutation_run:
	input:
		bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz",
		bed_index = "processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz.tbi",
		covariates = "processed/{study}/qtltools/input/{annot_type}/{condition}.covariates_prop.txt",
		vcf = config["qtl_vcf"]
	output:
		temp("processed/{study}/qtltools/output/{annot_type}/batches/{condition}.permutation.batch.{batch}.{n_batches}.txt")
	params:
		chunk = "{batch} {n_batches}"
	threads: 1
	resources:
		mem = 5000
	shell:
		"QTLtools_1.1_Ubuntu12.04_x86_64 cis --vcf {input.vcf} --bed {input.bed} --cov {input.covariates} --chunk {params.chunk} --out {output} --window {config[cis_window]} --permute 10000 --grp-best"


#Merge all batches from QTLtools
rule merge_permutation_batches:
	input:
		expand("processed/{{study}}/qtltools/output/{{annot_type}}/batches/{{condition}}.permutation.batch.{batch}.{n_batches}.txt", 
			batch=[i for i in range(1, config["n_batches"] + 1)],
			n_batches = config["n_batches"])
	output:
		"processed/{study}/qtltools/output/{annot_type}/{condition}.permuted.txt.gz"
	resources:
		mem = 100
	threads: 1
	shell:
		"cat {input} | bgzip > {output}"


#Run QTLtools in nominal mode
rule nominal_run:
	input:
		bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz",
		bed_index = "processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz.tbi",
		covariates = "processed/{study}/qtltools/input/{annot_type}/{condition}.covariates_prop.txt",
		vcf = config["qtl_vcf"]
	output:
		temp("processed/{study}/qtltools/output/{annot_type}/nominal_batches/{condition}.nominal.batch.{batch}.{n_batches}.txt")
	params:
		chunk = "{batch} {n_batches}"
	threads: 1
	resources:
		mem = 5000
	shell:
		"QTLtools_1.1_Ubuntu12.04_x86_64 cis --vcf {input.vcf} --bed {input.bed} --cov {input.covariates} --chunk {params.chunk} --out {output} --window {config[nominal_cis_window]} --nominal 1"

#Merge all batches from QTLtools
rule merge_nominal_batches:
	input:
		expand("processed/{{study}}/qtltools/output/{{annot_type}}/nominal_batches/{{condition}}.nominal.batch.{batch}.{n_batches}.txt", 
			batch=[i for i in range(1, config["n_batches"] + 1)],
			n_batches = config["n_batches"])
	output:
		temp("processed/{study}/qtltools/output/{annot_type}/{condition}.nominal.txt.gz")
	resources:
		mem = 100
	threads: 1
	shell:
		"cat {input} | bgzip > {output}"

#Add SNP coordinates to QTLTools output file
rule sort_qtltools_output:
	input:
		"processed/{study}/qtltools/output/{annot_type}/{condition}.nominal.txt.gz"
	output:
		protected("processed/{study}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz")
	resources:
		mem = 1000
	threads: 2
	shell:
		"zcat {input} | awk -v OFS='\\t' '{{$1=$1; print $0}}' | sort -k9,9 -k10,10n -k11,11n | bgzip > {output}"

#Tabix-index QTLtools output files
rule index_qtltools_output:
	input:
		"processed/{study}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz"
	output:
		"processed/{study}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz.tbi"
	resources:
		mem = 1000
	threads: 1
	shell:
		"tabix -s9 -b10 -e11 -f {input}"

#Run eigenMT correction on the linear model
rule run_eigenMT_chr:
	input:
		qtl = "processed/{study}/eigenMT/input/{condition}.input.txt",
		gen = "processed/{study}/eigenMT/genotypes/{chromosome}.genotypes.txt",
		gen_pos = "processed/{study}/eigenMT/genotypes/{chromosome}.snp_positions.txt",
		phe_pos = "processed/{study}//eigenMT/input/gene_positions.txt",
	params:
		chr = "{chromosome}"
	output:
		"processed/{study}/eigenMT/output/{condition}.{chromosome}.txt"
	resources:
		mem = 4000
	threads: 1
	shell:
		"python ~/software/eigenMT/eigenMT.py --QTL {input.qtl} --GEN {input.gen} --OUT {output} "
		"--GENPOS {input.gen_pos} --PHEPOS {input.phe_pos} --cis_dist {config[cis_window]} --external --CHROM {params.chr}"