#Lift over imputed genotypes from GRCh37 to GRCh38
cut -f1 genotypes/file_list.txt | python ~/software/utils/submitJobs.py --MEM 5000 --jobname liftOverVCF --command "python ~/software/utils/vcf/liftoverVcfGenotypes.py --chrMapFwd macrophage-gxe-study/data/liftOver_genotypes/GRCh38ToHg38_chromosome_map.txt --chrMapRev macrophage-gxe-study/data/liftOver_genotypes/Hg38ToGRCh38_chromosome_map.txt --liftOver macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain --reference ../../annotations/hg38/hg38.fa --vcfSuffix .vcf.gz --indir /nfs/users/nfs_k/ka8/group-scratch/hipsci/genotypes/REL-2014-11 --outdir genotypes/GRCh38/imputed_20151005 --execute True"

#### Make VCF file with SNPs only for rasqual ####

#Filter VCF files for samples, MAF
cut -f1 genotypes/file_list.txt | python ~/software/utils/submitJobs.py --MEM 1000 --jobname filterVCF --command "python ~/software/utils/vcf/filterVcf.py  --sampleList macrophage-gxe-study/data/sample_lists/SL1344/SL1344_gt_list.txt --MAF 0.05 --indir genotypes/GRCh38/imputed_20151005/ --outdir genotypes/SL1344/imputed_20151005/ --execute True --vcfSuffix .GRCh38.sorted.vcf.gz"

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
bcftools concat -a chr10.vcf.gz chr11.vcf.gz chr12.vcf.gz chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr1.vcf.gz chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz chr2.vcf.gz chr3.vcf.gz chr4.vcf.gz chr5.vcf.gz chr6.vcf.gz chr7.vcf.gz chr8.vcf.gz chr9.vcf.gz chrX.vcf.gz > imputed.86_samples.vcf

#Sort VCF files
#Sort genotypes (some chr1 snps went to chr9, GATK complains)
bsub -G team170 -n1 -R "span[hosts=1] select[mem>1000] rusage[mem=1000]" -q normal -M 1000 -o sortvcf.%J.jobout "~/software/vcflib/bin/vcfsort imputed.86_samples.vcf > imputed.86_samples.sorted.vcf"

#Run through vcfuniq and get rid of multiallelic SNPs, because liftover might have led to duplicated entries
~/software/vcflib/bin/vcfuniq imputed.86_samples.sorted.vcf | bcftools norm -m+any - | bcftools view -Oz -m2 -M2 - > imputed.86_samples.sorted.filtered.vcf.gz

#Add unique names to unnamed SNPs and remove duplicate snps
#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' imputed.86_samples.sorted.filtered.vcf.gz | bcftools view -O z -e 'ID=@duplicate_snps.txt' - > imputed.86_samples.sorted.filtered.named.vcf.gz
bcftools annotate --set-id +'%CHROM\_%POS' imputed.86_samples.sorted.filtered.vcf.gz | bcftools view -O z -e 'ID=@duplicate_snps.txt' - > imputed.86_samples.sorted.filtered.named.vcf.gz &

##### SNPS ONLY VCF for RASQUAL
bcftools view -v snps -O z imputed.86_samples.sorted.filtered.named.vcf.gz > imputed.86_samples.snps_only.vcf.gz &

#Extract SNP coords from a vcf file for RASQUAL
zgrep -v "#" imputed.86_samples.sorted.filtered.named.vcf.gz | cut -f 1,2,3 > imputed.86_samples.snp_coords.txt &
zgrep -v "#" imputed.86_samples.snps_only.vcf.gz | cut -f 1,2,3 > imputed.86_samples.snps_only.snp_coords.txt &

#Extract REF and ALT alleles and allele counts
zgrep -v "#" imputed.86_samples.sorted.filtered.named.vcf.gz | cut -f 1,2,3,4,5 > imputed.86_samples.variant_information.txt &
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\n' imputed.86_samples.sorted.filtered.named.vcf.gz | bgzip > imputed.86_samples.variant_information.txt.gz

#Filter by INFO score
bcftools filter -i 'INFO[0] >= 0.7' -O z imputed.86_samples.sorted.filtered.named.vcf.gz > imputed.86_samples.sorted.filtered.named.INFO_07.vcf.gz
zgrep -v "#" imputed.86_samples.sorted.filtered.named.INFO_07.vcf.gz | cut -f 1,2,3 > imputed.86_samples.snp_coords.INFO_07.txt &

#Split VCF file into chromosomes
cat macrophage-gxe-study/data/sample_lists/chromosome_list.txt | python ~/software/utils/vcf/vcfSplitByChromosome.py --vcf genotypes/SL1344/imputed_20151005/imputed.86_samples.sorted.filtered.named.vcf.gz --outdir genotypes/SL1344/imputed_20151005/chromosomes/ 

#Convert vcfs to GDS
/software/R-3.1.2/bin/Rscript ~/software/utils/vcf/vcfToGds.R --vcf-directory chromosomes_IFNO_07 --chr-list 

