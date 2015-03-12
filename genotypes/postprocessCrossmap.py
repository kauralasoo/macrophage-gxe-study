import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description = "Postprocess the VCF file created by the CrossMap to make it valid again.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--vcf", help = "Path to the VCF file.")
args = parser.parse_args()

vcf_file = open(args.vcf)
contigs = dict()
for line in vcf_file:
	line = line.rstrip()
	if(line[0] == "#"):
		print(line)
		if(line[0:8] == "##contig"):
			contig = line.split("##contig=<ID=")[1].split(",length=")[0]
			contigs[contig] = 1
	else:
		fields = line.split("\t",1)
		fields[0] = "chr" + fields[0] #Add chr in front of the chromosome names
		if(fields[0] in contigs):
			print("\t".join(fields)) #Only keep SNPs that fall into contigs mentioned in the header
