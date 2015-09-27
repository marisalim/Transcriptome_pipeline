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
library(vegan)
library(dplyr)
library(tidyr)

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
jpeg("LRTplot.jpg", height=5, width=5, units="in", res=500)
plot(rownums, LRTdat$LRT, pch=20, cex=1.5, xlab="Gene ID", ylab="LRT")
abline(h=3.84, lwd=2, col="tomato")
abline(h=6.64, lty=2, lwd=2, col="deepskyblue2")
abline(h=10.83, lty=4, lwd=2, col="seagreen3")
dev.off()

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

alldata_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006', 'namespace_1003'), filters='ensembl_peptide_id', values=LRTdat$Gene, mart=ensembl)
head(alldata_go) #i don't think HBB is here because it doesn't have a zebra finch ensembl id
alldata_BP <- alldata_go[alldata_go$namespace_1003 == "biological_process",] # not all ensembl ids have BP
alldata_MF <- alldata_go[alldata_go$namespace_1003 == "molecular_function",] # not all ensembl ids have MF
alldata_CC <- alldata_go[alldata_go$namespace_1003 == "cellular_component",] # not all ensembl ids have MF

head(alldata_MF)
head(alldata_BP)
head(alldata_CC)

alldata_MF$value <- 1
alldata_BP$value <- 1
alldata_CC$value <- 1

dim(alldata_MF)
dim(alldata_BP)
dim(alldata_CC)

LRTdat$signif_grp_05 <- ifelse(LRTdat$LRT > 3.84, "signif", "nonsig")
LRTdat$signif_grp_01 <- ifelse(LRTdat$LRT > 6.64, "signif", "nonsig")
LRTdat$signif_grp_001 <- ifelse(LRTdat$LRT > 10.83, "signif", "nonsig")
head(LRTdat)

MFmerge <- merge(x=alldata_MF, y=LRTdat, by.x="ensembl_peptide_id", by.y="Gene")
dim(MFmerge)
head(MFmerge)
BPmerge <- merge(x=alldata_BP, y=LRTdat, by.x="ensembl_peptide_id", by.y="Gene")
dim(BPmerge)
head(BPmerge)
CCmerge <- merge(x=alldata_CC, y=LRTdat, by.x="ensembl_peptide_id", by.y="Gene")
dim(CCmerge)
head(CCmerge)

MFdat <- spread(data=MFmerge, key=name_1006, value=value)
MFdat[is.na(MFdat)] <- 0
BPdat <- spread(data=BPmerge, key=name_1006, value=value)
BPdat[is.na(BPdat)] <- 0
CCdat <- spread(data=CCmerge, key=name_1006, value=value)
CCdat[is.na(CCdat)] <- 0
head(CCdat)

# Non-metric multidimensional scaling
MF_NMDS = metaMDS(MFdat[12:ncol(MFdat)], k=2)
BP_NMDS = metaMDS(BPdat[12:ncol(BPdat)], k=2)
CC_NMDS = metaMDS(CCdat[12:ncol(CCdat)], k=2)
# test different distance metrics

stressplot(MF_NMDS)
stressplot(BP_NMDS)
stressplot(CC_NMDS)


jpeg("all_NMDS.jpg", height=10, width=10, units="in", res=500)
par(mfrow=c(3,3))
ordiplot(MF_NMDS,type="n", main="MF, p=0.05")
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)

ordiplot(MF_NMDS,type="n", main="MF, p=0.01")
ordihull(MF_NMDS, groups=MFdat$signif_grp_01, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(MF_NMDS, groups=MFdat$signif_grp_01, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)

ordiplot(MF_NMDS,type="n", main="MF, p=0.001")
ordihull(MF_NMDS, groups=MFdat$signif_grp_001, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(MF_NMDS, groups=MFdat$signif_grp_001, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)

ordiplot(BP_NMDS,type="n", main="BP, p=0.05")
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)

ordiplot(BP_NMDS,type="n", main="BP, p=0.01")
ordihull(BP_NMDS, groups=BPdat$signif_grp_01, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(BP_NMDS, groups=BPdat$signif_grp_01, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)

ordiplot(BP_NMDS,type="n", main="BP, p=0.001")
ordihull(BP_NMDS, groups=BPdat$signif_grp_001, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(BP_NMDS, groups=BPdat$signif_grp_001, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)

ordiplot(CC_NMDS,type="n", main="CC, p=0.05")
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)

ordiplot(CC_NMDS,type="n", main="CC, p=0.01")
ordihull(CC_NMDS, groups=CCdat$signif_grp_01, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(CC_NMDS, groups=CCdat$signif_grp_01, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)

ordiplot(CC_NMDS,type="n", main="CC, p=0.001")
ordihull(CC_NMDS, groups=CCdat$signif_grp_001, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(CC_NMDS, groups=CCdat$signif_grp_001, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)
dev.off()

jpeg("0.05_NMDS.jpg", height=10, width=15, units="in", res=500)
par(mfrow=c(1,3))
ordiplot(MF_NMDS,type="n", main="MF, p=0.05")
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)

ordiplot(BP_NMDS,type="n", main="BP, p=0.05")
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)

ordiplot(CC_NMDS,type="n", main="CC, p=0.05")
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)
dev.off()

# individual plots for p=0.05
jpeg("0.05MF_NMDS.jpg", height=10, width=10, units="in", res=500)
ordiplot(MF_NMDS,type="n", main="MF, p=0.05")
orditorp(MF_NMDS,display="species",cex=1.5,col="tomato",air=0.01)
orditorp(MF_NMDS,display="sites",cex=1,air=0.1)
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)
dev.off()

jpeg("0.05BP_NMDS.jpg", height=10, width=10, units="in", res=500)
ordiplot(BP_NMDS,type="n", main="BP, p=0.05")
orditorp(BP_NMDS,display="species",cex=1.5,col="tomato",air=0.01)
orditorp(BP_NMDS,display="sites",cex=1,air=0.1)
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)
dev.off()

jpeg("0.05CC_NMDS.jpg", height=10, width=10, units="in", res=500)
ordiplot(CC_NMDS,type="n", main="CC, p=0.05")
orditorp(CC_NMDS,display="species",cex=1.5,col="tomato",air=0.01)
orditorp(CC_NMDS,display="sites",cex=1,air=0.1)
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="green")
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="blue", lty=2)
dev.off()



# -------------------- GO terms for each signif level ---------------------------
# p0.05_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=p0.05$Gene, mart=ensembl)
# head(p0.05_go)
# dim(p0.05_go)
# p0.05_BP <- p0.05_go2[p0.05_go2$namespace_1003 == "biological_process",]
# head(p0.05_BP)
# dat_tab <- table(p0.05_BP$ensembl_peptide_id, p0.05_BP$name_1006) 
# dat_tab[1:5, 1:5] 
#   #first column is an unidentified gene - that's why the colname is blank
#   #out of 80, there are only 60 with known biological_process functions
# write.csv(dat_tab, "dat_tab.csv")
# 
# nonsig_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006', 'namespace_1003'), filters='ensembl_peptide_id', values=nonsig_genes$Gene, mart=ensembl)
# head(nonsig_go)
# dim(nonsig_go)
# nonsig_BP <- nonsig_go[nonsig_go$namespace_1003 == "biological_process",]
# nonsig_tab <- table(nonsig_BP$ensembl_peptide_id, nonsig_BP$name_1006)
# nonsig_tab[1:5, 1:5]

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












