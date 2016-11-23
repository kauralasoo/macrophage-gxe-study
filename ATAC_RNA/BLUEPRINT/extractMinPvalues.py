import os
import sys
import argparse
import fileinput
import subprocess
import gzip

parser = argparse.ArgumentParser(description = "Extract minimal p-values per gene/peak from the BLUEPRINT summary stats.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--summary", help = "Path to the summary stats file.")
args = parser.parse_args()

summary_file = gzip.open(args.summary)
min_pvalue_dict = dict()
for line in summary_file:
	line = line.rstrip()
	fields = line.split(" ")
	phenotype = fields[2]
	pvalue = fields[3]
	if phenotype not in min_pvalue_dict:
		min_pvalue_dict[phenotype] = fields
	else:
		old_fields = min_pvalue_dict[phenotype]
		if(float(pvalue) < float(old_fields[3])):
			min_pvalue_dict[phenotype] = fields

#Print dict
for line in min_pvalue_dict.itervalues():
	print("\t".join(line))
