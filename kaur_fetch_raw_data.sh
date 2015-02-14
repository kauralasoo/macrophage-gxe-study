#Fetch all file names from iRODS
python ~/software/utils/irodsGetSamplesInStudy.py --studyName "Genetics of gene expression in human macrophage response to Salmonella" |  cut -f1 -d "." | uniq > fastq/SL1344_samples.txt

#Match file names to sample names
python ~/software/utils/irodsFetchMeta.py --irodsList fastq/SL1344_samples.txt | sort -k1 > fastq/SL1344_names.txt 

#Download fastq files from iRODS
cut -f1 fastq/SL1344_samples.txt | head -n 100 |  python ~/software/utils/irodsToFastq.py --output fastq/SL1344/
cut -f1 fastq/SL1344_samples.txt | head -n 200 | tail -n 100 |  python ~/software/utils/irodsToFastq.py --output fastq/SL1344/
cut -f1 fastq/SL1344_samples.txt | tail -n 328 |  python ~/software/utils/irodsToFastq.py --output fastq/SL1344/

