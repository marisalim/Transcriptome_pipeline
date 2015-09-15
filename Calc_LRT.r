# calculate LRT from codeml results
# null = Model A1 (NSsites = 2, model = 2, fix_omega = 1)
# test = Model A (NSsites = 2, model = 2, fix_omega = 0)

# Marisa Lim (c)2015

MLwd = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files"
setwd(MLwd)

nullmod <- read.table("lnL_null.txt")
testmod <- read.table("lnL_pos.txt")
head(nullmod)
head(testmod)

colnames(nullmod) <- c("Gene", "lnL")
colnames(testmod) <- c("Gene", "lnL")

LRTdat <- data.frame("Gene_nullmod"=nullmod$Gene, "lnL_nullmod"=nullmod$lnL, 
                     "Gene_testmod"=testmod$Gene, "lnL_testmod"=testmod$lnL)
LRTdat$LRT <- 2*(LRTdat$lnL_testmod - LRTdat$lnL_nullmod)

write.csv(LRTdat, "Loglikelihood_codemlresults.csv")
