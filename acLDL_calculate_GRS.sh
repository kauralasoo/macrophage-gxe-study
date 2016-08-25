#Extract genotypes from the hg19 files
cut -f1 genotypes/file_list.txt  | python ~/software/utils/submitJobs.py --MEM 1000 --jobname filterVCF --command "python ~/software/utils/vcf/filterVcf.py  --sampleList macrophage-gxe-study/data/sample_lists/acLDL/acLDL_gt_list.txt --indir /nfs/users/nfs_k/ka8/group-scratch/hipsci/genotypes/REL-2014-11 --outdir genotypes/acLDL/imputed_20151005/hg19 --execute True --vcfSuffix .vcf.gz"

#Rename files
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr10.filtered.vcf.gz chr10.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr11.filtered.vcf.gz chr11.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr12.filtered.vcf.gz chr12.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr13.filtered.vcf.gz chr13.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr14.filtered.vcf.gz chr14.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr15.filtered.vcf.gz chr15.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr16.filtered.vcf.gz chr16.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr17.filtered.vcf.gz chr17.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr18.filtered.vcf.gz chr18.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr19.filtered.vcf.gz chr19.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr1.filtered.vcf.gz chr1.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr20.filtered.vcf.gz chr20.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr21.filtered.vcf.gz chr21.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr22.filtered.vcf.gz chr22.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr2.filtered.vcf.gz chr2.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr3.filtered.vcf.gz chr3.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr4.filtered.vcf.gz chr4.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr5.filtered.vcf.gz chr5.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr6.filtered.vcf.gz chr6.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr7.filtered.vcf.gz chr7.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr8.filtered.vcf.gz chr8.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr9.filtered.vcf.gz chr9.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chrX.filtered.vcf.gz chrX.vcf.gz

#Index the vcf files
bcftools index chr10.vcf.gz
bcftools index chr11.vcf.gz
bcftools index chr12.vcf.gz
bcftools index chr13.vcf.gz
bcftools index chr14.vcf.gz
bcftools index chr15.vcf.gz
bcftools index chr16.vcf.gz
bcftools index chr17.vcf.gz
bcftools index chr18.vcf.gz
bcftools index chr19.vcf.gz
bcftools index chr1.vcf.gz
bcftools index chr20.vcf.gz
bcftools index chr21.vcf.gz
bcftools index chr22.vcf.gz
bcftools index chr2.vcf.gz
bcftools index chr3.vcf.gz
bcftools index chr4.vcf.gz
bcftools index chr5.vcf.gz
bcftools index chr6.vcf.gz
bcftools index chr7.vcf.gz
bcftools index chr8.vcf.gz
bcftools index chr9.vcf.gz
bcftools index chrX.vcf.gz

#Merge chromosomes
bcftools concat -a chr10.vcf.gz chr11.vcf.gz chr12.vcf.gz chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr1.vcf.gz chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz chr2.vcf.gz chr3.vcf.gz chr4.vcf.gz chr5.vcf.gz chr6.vcf.gz chr7.vcf.gz chr8.vcf.gz chr9.vcf.gz chrX.vcf.gz > imputed.70_samples.hg19.vcf

#Convert vcf into PLINK format
bsub -G team170 -n1 -R "span[hosts=1] select[mem>4000] rusage[mem=4000]" -q normal -M 4000 -o plink_convert.%J.jobout "plink --vcf imputed.70_samples.hg19.vcf --make-bed --out imputed.70_samples.hg19"

#Use PLINK to compute GRS
plink --bed imputed.70_samples.hg19.bed --score 79k_grs_weights.txt --extract 49k_grs_subset.txt --out grs_score.txt

