
#Map QTLs
snakemake --cluster ../../../../macrophage-trQTLs/scripts/snakemake_submit.py -np -s map_QTLs.snakefile processed/Fairfax/out.txt --jobs 1200 --configfile fairfax_config.yaml
