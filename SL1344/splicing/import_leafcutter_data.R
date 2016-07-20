library("readr")

#Import intron and cluster counts from leafcutter
intron_counts = readr::read_delim("results/SL1344/leafcutter/leafcutter.intron_counts.txt.gz", 
                                  delim = " ", col_names = TRUE)
cluster_counts = readr::read_delim("results/SL1344/leafcutter/leafcutter.cluster_counts.txt.gz", 
                                  delim = " ", col_names = TRUE)
