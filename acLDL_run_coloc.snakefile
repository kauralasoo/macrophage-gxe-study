configfile: "acLDL/acLDL_config.yaml"

rule run_all:
	input:
		expand("processed/acLDL/coloc/{gwas}.{phenotype}.{coloc_window}.txt", gwas = config["gwas_traits"], phenotype = config["coloc_phenotypes"], coloc_window = config["coloc_window"])
	output:
		"processed/acLDL/coloc_out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"

#Run coloc accross inflammatory traits
rule run_coloc:
	input:
		"processed/acLDL/fastqtl_output/{phenotype}/sorted/Ctrl.nominal.sorted.txt.gz"
	output:
		"processed/acLDL/coloc/{gwas}.{phenotype}.{coloc_window}.txt"
	params:
		outdir = "processed/acLDL/coloc",
		phenotype = "{phenotype}",
		gwas = "{gwas}",
		coloc_window = "{coloc_window}"
	resources:
		mem = 6000
	threads: 1
	shell:
		"/software/R-3.3.0/bin/Rscript acLDL/GWAS_overlaps/GWAS_run_coloc.R --phenotype {wildcards.phenotype} --window {wildcards.coloc_window} "
		"--gwas {wildcards.gwas} --dir {config[gwas_dir]} --outdir {params.outdir}"