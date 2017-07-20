#Process genotype data
snakemake --cluster ../../../../macrophage-trQTLs/scripts/snakemake_submit.py -np -s process_genotypes.snakefile processed/Fairfax/merged_genotypes/fairfax_genotypes.variant_information.txt.gz --jobs 100

#Convert the VCF file into an R matrix
bsub -G team170 -n1 -R "span[hosts=1] select[mem>35000] rusage[mem=35000]" -q normal -M 35000 -o importVCF.%J.jobout "/software/R-3.4.0/bin/Rscript macrophage-gxe-study/Fairfax/munge/importVCF.R"

#Map QTLs
snakemake --cluster ../../../../macrophage-trQTLs/scripts/snakemake_submit.py -np -s map_QTLs.snakefile processed/Fairfax/out.txt --jobs 1200 --configfile fairfax_config.yaml

#Run coloc against QTLs
snakemake --cluster ../macrophage-trQTLs/scripts/snakemake_submit.py -np -s macrophage-gxe-study/fairfax_run_coloc.snakefile processed/Fairfax/coloc_out.txt --jobs 100 --configfile fairfax_config.yaml

