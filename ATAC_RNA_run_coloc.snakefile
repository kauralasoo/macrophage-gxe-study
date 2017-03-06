configfile: "ATAC_RNA/ATAC_RNA_config.yaml"

rule run_all:
	input:
		expand("processed/ATAC_RNA/coloc/{gwas}.{phenotype}.{coloc_window}.txt", gwas = config["gwas_traits"], phenotype = config["coloc_phenotypes"], coloc_window = config["coloc_window"])
		expand("processed/ATAC_RNA/eQTL_caQTL_coloc/eQTL_caQTL_coloc.{condition}.{chromosome}.txt", condition = config["conditions"], chromosome = config["chromosomes"])
	output:
		"processed/ATAC_RNA/coloc_out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"

#Run coloc accross inflammatory traits
rule run_coloc:
	output:
		"processed/ATAC_RNA/coloc/{gwas}.{phenotype}.{coloc_window}.txt"
	params:
		outdir = "processed/ATAC_RNA/coloc",
		phenotype = "{phenotype}",
		gwas = "{gwas}",
		coloc_window = "{coloc_window}"
	resources:
		mem = 6000
	threads: 1
	shell:
		"/software/R-3.3.0/bin/Rscript ATAC_RNA/GWAS_overlaps/GWAS_run_coloc.R --phenotype {wildcards.phenotype} --window {wildcards.coloc_window} "
		"--gwas {wildcards.gwas} --dir {config[gwas_dir]} --outdir {params.outdir}"

#Run coloc accross inflammatory traits
rule run_caQTL_coloc:
	output:
		"processed/ATAC_RNA/eQTL_caQTL_coloc/eQTL_caQTL_coloc.{condition}.{chromosome}.txt"
	params:
		outdir = "processed/ATAC_RNA/eQTL_caQTL_coloc"
	resources:
		mem = 6000
	threads: 1
	shell:
		"/software/R-3.3.0/bin/Rscript ATAC_RNA/caQTLs_to_eQTLs/caQTL_eQTL_coloc.R --chr {wildcards.chromosome} "
		"--stimulation {wildcards.condition} --outdir {params.outdir}"