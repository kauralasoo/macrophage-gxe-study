#Extract genotypes from the cis regions of the three surface markers
bcftools view -r 1:161041759-162131963 -S macrophage-gxe-study/data/sample_lists/flow_cytometry_gt_list.txt genotypes/GRCh38/imputed_20151005/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr1.GRCh38.sorted.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - | bcftools norm -m+both - | bcftools view -m2 -M2 - > flow/genotypes/CD16_cis.vcf

bcftools view -r 5:140131728-141133701 -S macrophage-gxe-study/data/sample_lists/flow_cytometry_gt_list.txt genotypes/GRCh38/imputed_20151005/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr5.GRCh38.sorted.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - | bcftools norm -m+both - | bcftools view -m2 -M2 - > flow/genotypes/CD14_cis.vcf

bcftools view -r 10:17309344-18309344 -S macrophage-gxe-study/data/sample_lists/flow_cytometry_gt_list.txt genotypes/GRCh38/imputed_20151005/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr10.GRCh38.sorted.vcf.gz | bcftools filter -i 'MAF[0] >= 0.05' - | bcftools norm -m+both - | bcftools view -m2 -M2 - > flow/genotypes/CD206_cis.vcf

#Concat genotypes together
bcftools concat flow/genotypes/CD14_cis.vcf flow/genotypes/CD16_cis.vcf flow/genotypes/CD206_cis.vcf > flow/genotypes/flow_cis_regions.vcf 

#Map QTLs for surface markers
bgzip results/flow/fastqtl/mean_intensity.txt && tabix -f -p bed results/flow/fastqtl/mean_intensity.txt.gz
bgzip results/flow/fastqtl/random_intensity.txt && tabix -f -p bed results/flow/fastqtl/random_intensity.txt.gz

#Run fastqtl
cat results/flow/fastqtl/all_chunk_table.txt | python ~/software/utils/fastqtl/runFastQTL.py --vcf flow/genotypes/flow_cis_regions.filtered.vcf.gz --bed results/flow/fastqtl/mean_intensity.txt.gz --W 200000 --out results/flow/fastqtl/mean_out --execute True --permute 10000

cat results/flow/fastqtl/all_chunk_table.txt | python ~/software/utils/fastqtl/runFastQTL.py --vcf flow/genotypes/flow_cis_regions.filtered.vcf.gz --bed results/flow/fastqtl/random_intensity.txt.gz --W 200000 --out results/flow/fastqtl/random_out --execute True --permute 10000