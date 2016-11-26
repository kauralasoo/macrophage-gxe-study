#Sort and index summary stats
cat macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.txt | python ~/software/utils/tabixGWASSummaryStats.py --indir databases/GWAS/summary

#Convert summary stats into format suitable for GARFIELD
cat macrophage-gxe-study/data/gwas_catalog/GWAS_summary_stat_list.txt | python ~/software/utils/GWAS/extractGarfieldColumns.py --indir databases/GWAS/summary --outdir databases/garfield-data/pval/

#Extract variant coords from the GARFIELD annotation files
cut -f1 -d" " chr1 > variant_coords/chr1
cut -f1 -d" " chr10 > variant_coords/chr10
cut -f1 -d" " chr11 > variant_coords/chr11
cut -f1 -d" " chr12 > variant_coords/chr12
cut -f1 -d" " chr13 > variant_coords/chr13
cut -f1 -d" " chr14 > variant_coords/chr14
cut -f1 -d" " chr15 > variant_coords/chr15
cut -f1 -d" " chr16 > variant_coords/chr16
cut -f1 -d" " chr17 > variant_coords/chr17
cut -f1 -d" " chr18 > variant_coords/chr18
cut -f1 -d" " chr19 > variant_coords/chr19
cut -f1 -d" " chr2 > variant_coords/chr2
cut -f1 -d" " chr20 > variant_coords/chr20
cut -f1 -d" " chr21 > variant_coords/chr21
cut -f1 -d" " chr22 > variant_coords/chr22
cut -f1 -d" " chr3 > variant_coords/chr3
cut -f1 -d" " chr4 > variant_coords/chr4
cut -f1 -d" " chr5 > variant_coords/chr5
cut -f1 -d" " chr6 > variant_coords/chr6
cut -f1 -d" " chr7 > variant_coords/chr7
cut -f1 -d" " chr8 > variant_coords/chr8
cut -f1 -d" " chr9 > variant_coords/chr9
cut -f1 -d" " chrX > variant_coords/chrX