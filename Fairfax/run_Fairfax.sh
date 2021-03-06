#Process genotype data
snakemake --cluster ../../../../macrophage-trQTLs/scripts/snakemake_submit.py -np -s process_genotypes.snakefile processed/Fairfax/out.txt --jobs 100

#Convert the VCF file into an R matrix
bsub -G team170 -n1 -R "span[hosts=1] select[mem>35000] rusage[mem=35000]" -q normal -M 35000 -o importVCF.%J.jobout "/software/R-3.4.0/bin/Rscript macrophage-gxe-study/Fairfax/munge/importVCF.R"

#Map QTLs
snakemake --cluster ../../../../macrophage-trQTLs/scripts/snakemake_submit.py -np -s map_QTLs.snakefile processed/Fairfax/out.txt --jobs 1200 --configfile fairfax_config.yaml

#Run coloc against QTLs
snakemake --cluster ../macrophage-trQTLs/scripts/snakemake_submit.py -np -s macrophage-gxe-study/fairfax_run_coloc.snakefile processed/Fairfax/coloc_out.txt --jobs 100 --configfile macrophage-gxe-study/Fairfax/QTLs/fairfax_config.yaml

#Permute genotypes for one permutation run
bcftools query -l fairfax_genotypes.sorted.filtered.vcf.gz | sort -R > permuted_individuals.txt
bcftools reheader -s permuted_individuals.txt fairfax_genotypes.sorted.filtered.vcf.gz > fairfax_genotypes.sorted.filtered.permuted.vcf.gz
tabix -p vcf fairfax_genotypes.sorted.filtered.permuted.vcf.gz

#Map QTLs again using permuted genotypes
snakemake --cluster ../../../../macrophage-trQTLs/scripts/snakemake_submit.py -np -s map_QTLs.snakefile processed/Fairfax_permuted/out.txt --jobs 1200 --configfile fairfax_config_permuted.yaml
