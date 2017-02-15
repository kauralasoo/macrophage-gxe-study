configfile: "acLDL/acLDL_config.yaml"


rule map_qtls:
	input:
		expand("processed/acLDL/fastqtl_output/ensembl_87/{condition}.permuted.txt.gz", condition = config["conditions"])
	output:
		"processed/acLDL/fastqtl_output/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"

#Compres and index input bed file
rule compress_bed:
	input:
		bed = "processed/acLDL/fastqtl_splicing/ensembl_87/{condition}.norm_prop.txt"
	output:
		bed = "processed/acLDL/fastqtl_splicing/ensembl_87/{condition}.norm_prop.txt.gz",
		bed_index = "processed/acLDL/fastqtl_splicing/ensembl_87/{condition}.norm_prop.txt.gz.tbi"
	threads: 1
	resources:
		mem = 100
	shell:
		"bgzip {input.bed} && tabix -p bed {output.bed}"

#Run QTLtools in permutation mode
rule permutation_run:
	input:
		bed = "processed/acLDL/fastqtl_splicing/ensembl_87/{condition}.norm_prop.txt.gz",
		bed_index = "processed/acLDL/fastqtl_splicing/ensembl_87/AcLDL.norm_prop.txt.gz.tbi",
		covariates = "processed/acLDL/fastqtl_splicing/ensembl_87/{condition}.covariates_prop.txt",
		vcf = config["qtl_vcf"]
	output:
		temp("processed/acLDL/fastqtl_output/ensembl_87/batches/{condition}.permutation.batch.{batch}.{n_batches}.txt")
	params:
		chunk = "{batch} {n_batches}"
	threads: 1
	resources:
		mem = 5000
	shell:
		"QTLtools_1.1_Ubuntu12.04_x86_64 cis --vcf {input.vcf} --bed {input.bed} --cov {input.covariates} --chunk {params.chunk} --out {output} --window 100000 --permute 10000 --grp-best"


#Merge all batches from QTLtools
rule merge_qtl_batches:
	input:
		expand("processed/acLDL/fastqtl_output/ensembl_87/batches/{{condition}}.permutation.batch.{batch}.{n_batches}.txt", 
			batch=[i for i in range(1, config["n_batches"] + 1)],
			n_batches = config["n_batches"])
	output:
		"processed/acLDL/fastqtl_output/ensembl_87/{condition}.permuted.txt.gz"
	resources:
		mem = 100
	threads: 1
	shell:
		"cat {input} | bgzip > {output}"
