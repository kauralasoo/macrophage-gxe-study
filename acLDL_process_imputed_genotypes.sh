#Filter VCF files for samples, MAF > 0.05
cut -f1 genotypes/file_list.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname filterVCF --command "python ~/software/utils/vcf/filterVcf.py --sampleList macrophage-gxe-study/data/sample_lists/acLDL/acLDL_gt_list.txt --MAF 0.05 --indir genotypes/GRCh38/imputed_20151005/ --outdir genotypes/acLDL/imputed_20151005/ --execute True --vcfSuffix .GRCh38.sorted.vcf.gz"

#Rename files
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr10.filtered.GRCh38.sorted.vcf.gz chr10.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr11.filtered.GRCh38.sorted.vcf.gz chr11.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr12.filtered.GRCh38.sorted.vcf.gz chr12.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr13.filtered.GRCh38.sorted.vcf.gz chr13.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr14.filtered.GRCh38.sorted.vcf.gz chr14.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr15.filtered.GRCh38.sorted.vcf.gz chr15.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr16.filtered.GRCh38.sorted.vcf.gz chr16.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr17.filtered.GRCh38.sorted.vcf.gz chr17.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr18.filtered.GRCh38.sorted.vcf.gz chr18.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr19.filtered.GRCh38.sorted.vcf.gz chr19.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr1.filtered.GRCh38.sorted.vcf.gz chr1.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr20.filtered.GRCh38.sorted.vcf.gz chr20.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr21.filtered.GRCh38.sorted.vcf.gz chr21.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr22.filtered.GRCh38.sorted.vcf.gz chr22.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr2.filtered.GRCh38.sorted.vcf.gz chr2.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr3.filtered.GRCh38.sorted.vcf.gz chr3.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr4.filtered.GRCh38.sorted.vcf.gz chr4.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr5.filtered.GRCh38.sorted.vcf.gz chr5.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr6.filtered.GRCh38.sorted.vcf.gz chr6.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr7.filtered.GRCh38.sorted.vcf.gz chr7.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr8.filtered.GRCh38.sorted.vcf.gz chr8.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chr9.filtered.GRCh38.sorted.vcf.gz chr9.vcf.gz
mv hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.INFO_0.4_filtered.20151005.genotypes.chrX.filtered.GRCh38.sorted.vcf.gz chrX.vcf.gz


#Index VCF files
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

#Merge VCF files
bcftools concat -a chr10.vcf.gz chr11.vcf.gz chr12.vcf.gz chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr1.vcf.gz chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz chr2.vcf.gz chr3.vcf.gz chr4.vcf.gz chr5.vcf.gz chr6.vcf.gz chr7.vcf.gz chr8.vcf.gz chr9.vcf.gz chrX.vcf.gz > imputed.45_samples.vcf

#Sort VCF files
#Sort genotypes (some chr1 snps went to chr9, GATK complains)
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o sortvcf.%J.jobout "~/software/vcflib/bin/vcfsort imputed.45_samples.vcf > imputed.45_samples.sorted.vcf"

#Run through vcfuniq and get rid of multiallelic SNPs, because liftover might have led to duplicated entries
~/software/vcflib/bin/vcfuniq imputed.45_samples.sorted.vcf | bcftools norm -m+any - | bcftools view -m2 -M2 - > imputed.45_samples.sorted.filtered.vcf

#Filter by INFO score
bcftools filter -i 'INFO[0] >= 0.8' -O z imputed.45_samples.sorted.filtered.vcf > imputed.45_samples.snps_indels.INFO_08.vcf.gz 
