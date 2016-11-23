#Extract minimal p-values
python ../../macrophage-gxe-study/ATAC_RNA/BLUEPRINT/extractMinPvalues.py --summary mono_K4ME1_log2rpm_peer_10_all_summary.txt.20052016.gz | gzip > mono_K4ME1_min_pvalues.txt.gz 
python ../../macrophage-gxe-study/ATAC_RNA/BLUEPRINT/extractMinPvalues.py --summary mono_K27AC_log2rpm_peer_10_all_summary.txt.20052016.gz | gzip > mono_K27AC_min_pvalues.txt.gz
python ../../macrophage-gxe-study/ATAC_RNA/BLUEPRINT/extractMinPvalues.py --summary mono_gene_nor_combat_peer_10_all_summary.txt.20052016.gz | gzip > gene_min_pvalues.txt.gz

