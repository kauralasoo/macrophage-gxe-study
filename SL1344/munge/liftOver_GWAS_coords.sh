#IBD GWAS hits (trans-ethnic meta-analysis)
awk -v OFS='\t' '{print "chr"$2,$3,$3,$1,$4,$5,$7,$10,$11,$12}' IBD_trans_ethnic_association_summ_stats_b37.txt | tail -n+2  > IBD_trans_ethnic_bed.txt

CrossMap.py bed ../../../../macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain IBD_trans_ethnic_bed.txt IBD_trans_ethnic_bed.Hg38.txt
awk '{sub(/chr/,"")}1' IBD_trans_ethnic_bed.Hg38.txt > IBD_trans_ethnic_bed.GRCh38.txt

sort -k1,1 -k2,2n IBD_trans_ethnic_bed.GRCh38.txt > IBD_trans_ethnic_bed.GRCh38.sorted.txt

#Convert to GWAS format for IGV
echo -e "CHR\tBP\tSNP\tP" > IBD_trans_ethnic.gwas && cut -f 1,2,4,10 IBD_trans_ethnic_bed.GRCh38.sorted.txt >> IBD_trans_ethnic.gwas


#Full GWAS p-values
#### IBD ###
#GWAS to BED
awk -v OFS='\t' '{print "chr"$1,$3,$3,$2,$4,$5,$8,$9,$10,$11}' EUR.IBD.gwas.assoc | tail -n+2 > EUR.IBD.gwas.bed.txt

#CrossMap
CrossMap.py bed ../../../macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain EUR.IBD.gwas.bed.txt EUR.IBD.gwas.bed.Hg38.txt
awk '{sub(/chr/,"")}1' EUR.IBD.gwas.bed.Hg38.txt > EUR.IBD.gwas.bed.GRCh38.txt

#Sort and index
sort -k1,1 -k2,2n EUR.IBD.gwas.bed.GRCh38.txt | bgzip > EUR.IBD.gwas.bed.GRCh38.sorted.txt.gz
tabix -p bed EUR.IBD.gwas.bed.GRCh38.sorted.txt.gz

#### CD ####
#GWAS to BED
awk -v OFS='\t' '{print "chr"$1,$3,$3,$2,$4,$5,$8,$9,$10,$11}' EUR.CD.gwas.assoc | tail -n+2 > EUR.CD.gwas.bed.txt

#CrossMap
CrossMap.py bed ../../../macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain EUR.CD.gwas.bed.txt EUR.CD.gwas.bed.Hg38.txt
awk '{sub(/chr/,"")}1' EUR.CD.gwas.bed.Hg38.txt > EUR.CD.gwas.bed.GRCh38.txt

#Sort and index
sort -k1,1 -k2,2n EUR.CD.gwas.bed.GRCh38.txt | bgzip > EUR.CD.gwas.bed.GRCh38.sorted.txt.gz
tabix -p bed EUR.CD.gwas.bed.GRCh38.sorted.txt.gz


#### UC ####
awk -v OFS='\t' '{print "chr"$1,$3,$3,$2,$4,$5,$8,$9,$10,$11}' EUR.UC.gwas.assoc | tail -n+2 > EUR.UC.gwas.bed.txt

#CrossMap
CrossMap.py bed ../../../macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain EUR.UC.gwas.bed.txt EUR.UC.gwas.bed.Hg38.txt
awk '{sub(/chr/,"")}1' EUR.UC.gwas.bed.Hg38.txt > EUR.UC.gwas.bed.GRCh38.txt

#Sort and index
sort -k1,1 -k2,2n EUR.UC.gwas.bed.GRCh38.txt | bgzip > EUR.UC.gwas.bed.GRCh38.sorted.txt.gz
tabix -p bed EUR.UC.gwas.bed.GRCh38.sorted.txt.gz


#Lift over IGAP GWAS coordinates
awk -v OFS="\t" '{print "chr"$1,$2,$2,$3,$4,$5,"NA",$6,$7,$8}' IGAP_stage_1.txt | tail -n+2 > IGAP_stage_1.bed
CrossMap.py bed ../../../macrophage-gxe-study/data/liftOver_genotypes/hg19ToHg38.over.chain IGAP_stage_1.bed IGAP_stage_1.Hg38.bed
awk '{sub(/chr/,"")}1' IGAP_stage_1.Hg38.bed > IGAP_stage_1.GRCh38.bed
sort -k1,1 -k2,2n IGAP_stage_1.GRCh38.bed | bgzip > IGAP_stage_1.GRCh38.sorted.bed.gz
tabix -p bed IGAP_stage_1.GRCh38.sorted.bed.gz




#Extract a region
echo -e "CHR\tBP\tSNP\tP" > IBD.RGS14.gwas && tabix EUR.IBD.gwas.bed.GRCh38.sorted.txt.gz 5:177,000,000-177,500,000 | cut -f 1,2,4,10 >> IBD.RGS14.gwas
echo -e "CHR\tBP\tSNP\tP" > CD.RGS14.gwas && tabix EUR.CD.gwas.bed.GRCh38.sorted.txt.gz 5:177,000,000-177,500,000 | cut -f 1,2,4,10 >> CD.RGS14.gwas
echo -e "CHR\tBP\tSNP\tP" > UC.RGS14.gwas && tabix EUR.UC.gwas.bed.GRCh38.sorted.txt.gz 5:177,000,000-177,500,000 | cut -f 1,2,4,10 >> UC.RGS14.gwas


echo -e "CHR\tBP\tSNP\tP" > IBD.UBQLN4.gwas && tabix EUR.IBD.gwas.bed.GRCh38.sorted.txt.gz 1:155536176-156553003 | cut -f 1,2,4,10 >> IBD.UBQLN4.gwas
echo -e "CHR\tBP\tSNP\tP" > IBD.chr1.gwas && tabix EUR.IBD.gwas.bed.GRCh38.sorted.txt.gz 1 | cut -f 1,2,4,10 >> IBD.chr1.gwas
