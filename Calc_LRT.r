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

dim(nullmod)
dim(testmod)

# get gene name
replacebreak <- gsub(pattern="/", replacement=" ", x=nullmod$Gene)
splitbreak <- strsplit(replacebreak, split=" ")
genename1 <- sapply(splitbreak, function(x){
  paste(x[[3]])
})
head(genename1)
replacebreak2 <- gsub(pattern="_", replacement=" ", x=genename1)
splitbreak2 <- strsplit(replacebreak2, split=" ")
genename <- sapply(splitbreak2, function(x){
  paste(x[[1]])
})
head(genename)

LRTdat <- data.frame("Gene"=genename, "Gene_nullmod"=nullmod$Gene, "lnL_nullmod"=nullmod$lnL, 
                     "Gene_testmod"=testmod$Gene, "lnL_testmod"=testmod$lnL)
LRTdat$LRT <- 2*(LRTdat$lnL_testmod - LRTdat$lnL_nullmod)
head(LRTdat)

# Which columns have significant LRT? (df = 1)
dim(LRTdat[LRTdat$LRT > 3.84,]) # p=0.05
dim(LRTdat[LRTdat$LRT > 6.64,]) # p=0.01
dim(LRTdat[LRTdat$LRT > 10.83,]) # p=0.01


# add column with gene name and/or GO terms...KEGG? 




write.csv(LRTdat, "Loglikelihood_codemlresults.csv")
