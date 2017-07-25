
#Subset
bgzip processed/Fairfax/qtltools/input/shared_84/CD14.aFC.txt && tabix -p bed processed/Fairfax/qtltools/input/shared_84/CD14.aFC.txt.gz
bgzip processed/Fairfax/qtltools/input/shared_84/IFN.aFC.txt && tabix -p bed processed/Fairfax/qtltools/input/shared_84/IFN.aFC.txt.gz
bgzip processed/Fairfax/qtltools/input/shared_84/LPS2.aFC.txt && tabix -p bed processed/Fairfax/qtltools/input/shared_84/LPS2.aFC.txt.gz
bgzip processed/Fairfax/qtltools/input/shared_84/LPS24.aFC.txt && tabix -p bed processed/Fairfax/qtltools/input/shared_84/LPS24.aFC.txt.gz

#Run aFC calculation
python ~/software/aFC/aFC.py --vcf processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz --pheno processed/Fairfax/qtltools/input/shared_84/CD14.aFC.txt.gz --qtl processed/Fairfax/qtltools/input/shared_84/qtl_pairs.txt --log_xform 1 --log_base 2 --o processed/Fairfax/qtltools/input/shared_84/CD14.aFC_results.txt

python ~/software/aFC/aFC.py --vcf processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz --pheno processed/Fairfax/qtltools/input/shared_84/IFN.aFC.txt.gz --qtl processed/Fairfax/qtltools/input/shared_84/qtl_pairs.txt --log_xform 1 --log_base 2 --o processed/Fairfax/qtltools/input/shared_84/IFN.aFC_results.txt

python ~/software/aFC/aFC.py --vcf processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz --pheno processed/Fairfax/qtltools/input/shared_84/LPS2.aFC.txt.gz --qtl processed/Fairfax/qtltools/input/shared_84/qtl_pairs.txt --log_xform 1 --log_base 2 --o processed/Fairfax/qtltools/input/shared_84/LPS2.aFC_results.txt

python ~/software/aFC/aFC.py --vcf processed/Fairfax/merged_genotypes/fairfax_genotypes.sorted.filtered.vcf.gz --pheno processed/Fairfax/qtltools/input/shared_84/LPS24.aFC.txt.gz --qtl processed/Fairfax/qtltools/input/shared_84/qtl_pairs.txt --log_xform 1 --log_base 2 --o processed/Fairfax/qtltools/input/shared_84/LPS24.aFC_results.txt


#Full dataset
bgzip processed/Fairfax/qtltools/input/shared/CD14.aFC.txt && tabix -p bed processed/Fairfax/qtltools/input/shared/CD14.aFC.txt.gz
bgzip processed/Fairfax/qtltools/input/shared/IFN.aFC.txt && tabix -p bed processed/Fairfax/qtltools/input/shared/IFN.aFC.txt.gz
bgzip processed/Fairfax/qtltools/input/shared/LPS2.aFC.txt && tabix -p bed processed/Fairfax/qtltools/input/shared/LPS2.aFC.txt.gz
bgzip processed/Fairfax/qtltools/input/shared/LPS24.aFC.txt && tabix -p bed processed/Fairfax/qtltools/input/shared/LPS24.aFC.txt.gz

#Calculate aFC
