library("purrr")
library("devtools")
library("rtracklayer")
load_all("../reviseAnnotations/")

#Import revised GFF files
gff_list = list(upstream = rtracklayer::import.gff3("processed/acLDL/annotations/gff/reviseAnnotations_upstream.gff3"),
                contained = rtracklayer::import.gff3("processed/acLDL/annotations/gff/reviseAnnotations_contained.gff3"),
                downstream = rtracklayer::import.gff3("processed/acLDL/annotations/gff/reviseAnnotations_downstream.gff3"))

#Convert the GFF files into GRanges lists
granges_lists = purrr::map(gff_list, ~reviseAnnotations::revisedGffToGrangesList(.))
saveRDS(granges_lists, "results/reviseAnnotations/reviseAnnotations.GRangesList.rds")
