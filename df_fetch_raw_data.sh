#20/08/14 - Download raw data from IRODS
python ~/software/utils/fetch-irods.py --runid 13569 --laneids 5,6,7 --sampleids 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --dir raw
python ~/software/utils/fetch-irods.py --runid 13518 --laneids 8 --sampleids 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --dir raw
python ~/software/utils/renameFiles.py --samples macrophage-gxe-study/data/DF_rename_runs_patch1.txt --filedir raw --suffix .bam --execute True
