#Rename bams - 15/09/14
python ~/software/utils/renameFiles.py --samples macrophage-gxe-study/data/ka-samples-batch_1.txt --filedir bams_merged --suffix .bam --execute False

#Reheader bam files - 15/09/14
cut -f2 macrophage-gxe-study/data/ka-samples-batch_1.txt | python ~/software/utils/reheaderBams.py --indir bams_merged/ --outdir bams_reheader/ --new_header bams_merged/old_header.txt --execute True

#Get preliminary read counts using featureCounts - 15/09/14
bsub -G team170 -o FarmOut/merge_bams.%J.txt "featureCounts -p -a ../../annotations/GRCh37/Ensembl_74/Homo_sapiens.GRCh37.74.gff -o results/ka-batch1.counts bams_reheader/ffdp_A.bam bams_reheader/ougl_A.bam bams_reheader/oomz_A.bam bams_reheader/liun_A.bam bams_reheader/ffdp_B.bam bams_reheader/ougl_B.bam bams_reheader/oomz_B.bam bams_reheader/liun_B.bam bams_reheader/ffdp_C.bam bams_reheader/ougl_C.bam bams_reheader/oomz_C.bam bams_reheader/liun_C.bam bams_reheader/ffdp_D.bam bams_reheader/ougl_D.bam bams_reheader/oomz_D.bam bams_reheader/liun_D.bam"