# calculate LRT from codeml results
# null = Model A1 (NSsites = 2, model = 2, fix_omega = 1)
# test = Model A (NSsites = 2, model = 2, fix_omega = 0)

# Marisa Lim (c)2015

# set working directory
MLwd = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files"
setwd(MLwd)

# load libraries
library(biomaRt)
library("GO.db")
#library(KEGGREST)
library(goProfiles)

# Read in lnL values for each model
nullmod <- read.table("lnL_null.txt")
testmod <- read.table("lnL_pos.txt")
head(nullmod)
head(testmod)
# give column names
colnames(nullmod) <- c("Gene", "lnL")
colnames(testmod) <- c("Gene", "lnL")
# check dimensions
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

# merge the information for both models
LRTdat <- data.frame("Gene"=genename, "Gene_nullmod"=nullmod$Gene, "lnL_nullmod"=nullmod$lnL, 
                     "Gene_testmod"=testmod$Gene, "lnL_testmod"=testmod$lnL)
# calculate LRT
LRTdat$LRT <- 2*(LRTdat$lnL_testmod - LRTdat$lnL_nullmod)
head(LRTdat)

# plot LRT and significance cutoffs
rownums <- c(1:163)
plot(rownums, LRTdat$LRT, pch=20, cex=1.5, xlab="Gene ID", ylab="LRT")
abline(h=3.84, lwd=2, col="tomato")
abline(h=6.64, lty=2, lwd=2, col="deepskyblue2")
abline(h=10.83, lty=4, lwd=2, col="seagreen3")

# Which columns have significant LRT? (df = 1)
dim(LRTdat[LRTdat$LRT >= 3.84,]) # p=0.05
dim(LRTdat[LRTdat$LRT >= 6.64,]) # p=0.01
dim(LRTdat[LRTdat$LRT >= 10.83,]) # p=0.001

# new dataframes for each level of significance
p0.05 <- LRTdat[LRTdat$LRT >= 3.84,]
p0.01 <- LRTdat[LRTdat$LRT >= 6.64,]
p0.001 <- LRTdat[LRTdat$LRT >= 10.83,]

# read in list of genes not significant in codeml (at any level)
nonsig_genes <- LRTdat[LRTdat$LRT < 3.84,]

# search GO terms/info
ensembl = useMart("ensembl", dataset="tguttata_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

p0.05_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=p0.05$Gene, mart=ensembl)
head(p0.05_go)
dim(p0.05_go)
p0.05_go2 <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006', 'namespace_1003'), filters='ensembl_peptide_id', values=p0.05$Gene, mart=ensembl)

p0.05_go2[p0.05_go2$hgnc_symbol == "TPP1", ]

p0.05_BP <- p0.05_go2[p0.05_go2$namespace_1003 == "biological_process",]
head(p0.05_BP)
dat_tab <- table(p0.05_BP$hgnc_symbol, p0.05_BP$name_1006) 
write.table(dat_tab, "p0.05gotermstable.txt")
barplot(dat_tab)

nonsig_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006', 'namespace_1003'), filters='ensembl_peptide_id', values=nonsig_genes$Gene, mart=ensembl)
head(nonsig_go)
dim(nonsig_go)
nonsig_BP <- nonsig_go[nonsig_go$namespace_1003 == "biological_process",]
nonsig_tab <- table(nonsig_BP$hgnc_symbol, nonsig_BP$name_1006)
barplot(nonsig_tab)


# p0.01_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=p0.01$Gene, mart=ensembl)
# head(p0.01_go)
# dim(p0.01_go)
# p0.01_go2 <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'go_id'), filters='ensembl_peptide_id', values=p0.01$Gene, mart=ensembl)
# p0.01goterms <- data.frame("go_description"=Term(p0.01_go2$go_id))
# head(p0.01goterms)
# 
# p0.001_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=p0.001$Gene, mart=ensembl)
# head(p0.001_go)
# dim(p0.001_go)
# p0.001_go2 <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'go_id'), filters='ensembl_peptide_id', values=p0.001$Gene, mart=ensembl)
# p0.001goterms <- data.frame("go_description"=Term(p0.001_go2$go_id))
# head(p0.001goterms)

# goal: is there enrichment of certain go categories for biological process?
# how to separate molecular function vs. biological process vs. cellular component go terms??

#### Try goProfiles?












