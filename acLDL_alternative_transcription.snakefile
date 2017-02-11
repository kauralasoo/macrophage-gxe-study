configfile: "acLDL/acLDL_config.yaml"

#Quantify gene expression using full Ensembl annotations
rule ensembl_quant_salmon:
	input:
		fq1 = "processed/acLDL/fastq/{sample}.1.fastq.gz",
		fq2 = "processed/acLDL/fastq/{sample}.2.fastq.gz"
	output:
		"processed/acLDL/salmon/ensembl_87/{sample}/quant.sf"
	params:
		out_prefix = "processed/acLDL/salmon/ensembl_87/{sample}"
	resources:
		mem = 10000
	threads: 8	
	shell:
		"salmon --no-version-check quant --useVBOpt --seqBias --gcBias --libType ISR "
		"--index {config[salmon_ensembl_index]} -1 {input.fq1} -2 {input.fq2} -p {threads} "
		"--geneMap {config[ensembl_gtf]} -o {params.out_prefix}"

#Align reads to the reference genome using STAR
rule star_align:
	input:
		fq1 = "processed/acLDL/fastq/{sample}.1.fastq.gz",
		fq2 = "processed/acLDL/fastq/{sample}.2.fastq.gz"
	output:
		bam = "processed/acLDL/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	params:
		prefix = "processed/acLDL/STAR/{sample}/{sample}.",
		rg = "ID:1\tSM:{sample}"
	resources:
		mem = 42000
	threads: 8
	shell:
		"STAR --runThreadN {threads} --outSAMtype BAM SortedByCoordinate --outWigType bedGraph "
		"--outWigNorm None --outWigStrand Stranded --outSAMattrRGline \"{params.rg}\" "
		"--readFilesCommand zcat --genomeDir {config[star_index]} --limitBAMsortRAM 10000000000 "
		"--outFileNamePrefix {params.prefix} --readFilesIn {input.fq1} {input.fq2} "

#Index sorted bams
rule index_bams:
	input:
		"processed/acLDL/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	output:
		"processed/acLDL/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
	resources:
		mem = 50
	threads: 1
	shell:
		"samtools index {input}"

#Check genotype concordance between RNA-seq and VCF
rule check_genotype_concordance:
	input:
		"processed/acLDL/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
		"processed/acLDL/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
	output:
		"processed/acLDL/verifyBamID/{sample}.verifyBamID.bestSM"
	params:
		out_prefix = "processed/acLDL/verifyBamID/{sample}.verifyBamID"
	resources:
		mem = 1500
	threads: 1
	shell:
		"verifyBamID.1.1.2 --vcf {config[vcf_file]} --bam {input} --out {params.out_prefix} --best --ignoreRG"

#Convert bedgraph to bigwig
rule bedgraph_to_bigwig:
	input:
		"processed/acLDL/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	output:
		bw1 = "processed/acLDL/bigwig/{sample}.str1.bw",
		bw2 = "processed/acLDL/bigwig/{sample}.str2.bw"
	params:
		bg1 = "processed/acLDL/STAR/{sample}/{sample}.Signal.Unique.str1.out.bg",
		bg2 = "processed/acLDL/STAR/{sample}/{sample}.Signal.Unique.str2.out.bg",
		bg3 = "processed/acLDL/STAR/{sample}/{sample}.Signal.UniqueMultiple.str1.out.bg",
		bg4 = "processed/acLDL/STAR/{sample}/{sample}.Signal.UniqueMultiple.str2.out.bg"
	resources:
		mem = 1000
	threads: 1
	shell:
		"bedGraphToBigWig {params.bg1} {config[chromosome_lengths]} {output.bw1} && "
		"bedGraphToBigWig {params.bg2} {config[chromosome_lengths]} {output.bw2} && "
		"rm {params.bg1} && rm {params.bg2} && rm {params.bg3} && rm {params.bg4}"




