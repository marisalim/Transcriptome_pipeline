# for plotting the Bayesian Empirical Bayes Site class proportions
# (c)2015 Marisa Lim, Stony Brook University

MLwd = "C:/Users/mcwlim/Desktop/StonyBrook/GrahamLab/Dissertation idea materials/THESIS PROJECTS/Hummingbird_Genomics/Transcriptome pipeline files/Candidate gene files/CladeC_results"
setwd(MLwd)

HBB <- read.csv("HBB_cladeC_BEB.csv", header=T)
head(HBB)
SCP2 <- read.csv("SCP2_cladeC_BEB.csv", header=T)
head(SCP2)
