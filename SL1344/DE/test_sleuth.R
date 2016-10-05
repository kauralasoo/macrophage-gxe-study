library('seluth')
library('wasabi')
library("devtools")
library("cqn")
library("dplyr")
load_all("../seqUtils/")

#Import sailfish quants
sf_dirs = dir("results/PRJDB2508_quants/")
sf_dirs = paste0("results/PRJDB2508_quants/", sf_dirs)

#Convert quants format suitable for sleuth
prepare_fish_for_sleuth(sf_dirs)

# we have 4 conditions: wild type, each gene knocked out individually,	 	 
# and both kocked out together	 	 
genotypes <- c(	 	 
  rep('wildtype', 4),	 	 
  rep('ros1-3', 4),	 	 
  rep('dml2;3', 4),	 	 
  rep('ros1dml2;3', 4)	 	 
)	 	 

# describe the experiment in terms of genotype and growth condition	 	 
# we also summarise the gene knockouts and osmotic stress status	 	 
# or each sample	 	 
s2c <- data.frame(	 	 
  sample = sf_dirs,	 	 
  genotype = genotypes,	 	 
  growth = rep(c(':mock', ':aba', ':saline', ':dry'), 4),	 	 
  ros1_3 = rep(c(rep(':on', 4), rep(':off', 4)), 2),	 	 
  dml2_3 = c(rep(':on', 8), rep(':off', 8)),	 	 
  osmotic_stress = rep(c(':no', ':no', ':yes', ':yes'), 4),
  path = sf_dirs
)	 	 
s2c$path <- as.character(s2c$path)	 	 
s2c	 

# model the experimental question	 	 
model <- "~ osmotic_stress"	 	 
# load data and fit the model	 	 
so <- sleuth_prep(s2c, as.formula(model)) %>%
  sleuth_fit()
models(so)

# test 
so <- sleuth_wt(so, which_beta = 'osmotic_stress:yes')	 	 

# inspect the results	 	 
sleuth_live(so)

#Extract gene lists
results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
