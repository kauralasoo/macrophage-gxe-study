configfile: "acLDL/acLDL_config.yaml"

#Quantify gene expression using full Ensembl annotations
rule ensmembl_quant_salmon:
	input:
		fq1 = "processed/acLDL/fastq/{sample}.1.fastq.gz",
		fq2 = "processed/acLDL/fastq/{sample}.2.fastq.gz"
	output:
		"processed/acLDL/salmon/ensembl_87/{sample}"
	resources:
		mem = 10000
	threads: 8	
	shell:
		"salmon --no-version-check quant --useVBOpt --seqBias --gcBias --libType ISR "
		"--index {config[salmon_ensembl_index]} -1 {input.fq1} -2 {input.fq2} -p {threads} "
		"--geneMap {config[ensembl_gff3]} -o {output}"