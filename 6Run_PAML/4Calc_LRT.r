# calculate LRT from codeml results, grab gene ontology categories, NMDS
# null = Model A1 (NSsites = 2, model = 2, fix_omega = 1)
# test = Model A (NSsites = 2, model = 2, fix_omega = 0)

# Marisa Lim (c)2015

# load libraries
library(biomaRt)
library("GO.db")
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# set working directory
MLwd = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files"
setwd(MLwd)

# ----------------- Read in lnL values for each model -----------------
# nuclear dna files
nu_nullmod <- read.table("nu_lnL_null.txt")
nu_testmod <- read.table("nu_lnL_pos.txt")
head(nu_nullmod)
head(nu_testmod)
# give column names
colnames(nu_nullmod) <- c("Gene", "lnL")
colnames(nu_testmod) <- c("Gene", "lnL")
# check dimensions
dim(nu_nullmod)
dim(nu_testmod)
# get gene name
nu_replacebreak <- gsub(pattern="/", replacement=" ", x=nu_nullmod$Gene)
nu_splitbreak <- strsplit(nu_replacebreak, split=" ")
nu_genename1 <- sapply(nu_splitbreak, function(x){
  paste(x[[7]])
})
head(nu_genename1)
nu_replacebreak2 <- gsub(pattern="_", replacement=" ", x=nu_genename1)
nu_splitbreak2 <- strsplit(nu_replacebreak2, split=" ")
nu_genename <- sapply(nu_splitbreak2, function(x){
  paste(x[[1]])
})
head(nu_genename)

# mitochondrial dna files
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
  paste(x[[7]])
})
head(genename1)
replacebreak2 <- gsub(pattern="_", replacement=" ", x=genename1)
splitbreak2 <- strsplit(replacebreak2, split=" ")
genename <- sapply(splitbreak2, function(x){
  paste(x[[1]])
})
head(genename)

# merge the information for both models 
LRTdat_nu <- data.frame("Gene"=nu_genename, "Gene_nullmod"=nu_nullmod$Gene, "lnL_nullmod"=nu_nullmod$lnL, 
                          "Gene_testmod"=nu_testmod$Gene, "lnL_testmod"=nu_testmod$lnL)
head(LRTdat_nu)
LRTdat_mito <- data.frame("Gene"=genename, "Gene_nullmod"=nullmod$Gene, "lnL_nullmod"=nullmod$lnL, 
                     "Gene_testmod"=testmod$Gene, "lnL_testmod"=testmod$lnL)
head(LRTdat_mito)

# combine nuclear and mito datasets
LRTdat <- rbind(LRTdat_nu, LRTdat_mito)
dim(LRTdat)
head(LRTdat)

# ----------------- calculate LRT -----------------
# Want values > 0, ignore any < 0
LRTdat$LRT <- 2*(LRTdat$lnL_testmod - LRTdat$lnL_nullmod)
head(LRTdat)

# plot LRT and significance cutoffs
LRTdat$rownums <- c(1:nrow(LRTdat))

# add signif level identifiers
signif_lev <- c()
for(i in 1:nrow(LRTdat)){
  checkthis <- LRTdat$LRT[i]
  
  if(checkthis > 10.83){
    signif_lev[i] <- 'signif001'
  } else if(checkthis > 6.64){
    signif_lev[i] <- 'signif01'
  } else if(checkthis > 3.84){
    signif_lev[i] <- 'signif05'
  } else{
    signif_lev[i] <- 'nonsig'
  }
}
signif_lev

LRTdat$Significance_level <- signif_lev
names(LRTdat)

ggplot(LRTdat, aes(x=rownums, y=LRT, color=Significance_level)) + 
  geom_point(size=3, alpha=0.6) + 
  theme_bw() + xlab("Index") + ylab("Likelihood ratio test value")
ggsave("LRTplot.jpg", height=10, width=12, units="in", dpi=500)

LRTdat2 <- LRTdat[LRTdat$LRT > 3.84,]
dim(LRTdat2)
p0.05 <- LRTdat[LRTdat$LRT >= 3.84,]

# if usual host not working today...25April16, the above host is an archived
#ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="tguttata_gene_ensembl", host="dec2015.archive.ensembl.org")
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="tguttata_gene_ensembl", host='www.ensembl.org')            

p0.05_go_names <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=p0.05$Gene, mart=ensembl)
LRTdat2 <- merge(x=LRTdat2, y=p0.05_go_names, by.x='Gene', by.y='ensembl_peptide_id')
LRTdat2$hgnc_symbol[17] <- 'COXI'
LRTdat2$hgnc_symbol
# remove messed up LRT genes: Cactin, EEF2, RPS4X, MTX2, COX1
LRTdat2 <- LRTdat2[-c(1:2, 7, 11, 17),]

ggplot(LRTdat2, aes(x=rownums, y=LRT)) + 
  geom_point(size=5, color='grey') + 
  geom_label_repel(aes(label=hgnc_symbol, fill=factor(Significance_level)), 
                  nudge_x=1, nudge_y=1, fontface='bold', color='white',
                  box.padding=unit(0.25, 'lines')) +
  theme_bw() + xlab("Index") + ylab("Likelihood ratio test value") +
  theme(text=element_text(size=20), legend.position='bottom')
ggsave("LRTplot_signifgenes.jpg", height=10, width=12, units="in", dpi=500)

# Which columns have significant LRT? (df = 1)
p05dat <- LRTdat[LRTdat$LRT >= 3.84,] # p=0.05
p01dat <- LRTdat[LRTdat$LRT >= 6.64,] # p=0.01
p001dat <- LRTdat[LRTdat$LRT >= 10.83,] # p=0.001
dim(p05dat)
dim(p01dat)
dim(p001dat)

p0.05_go_info <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006','go_id', 'namespace_1003'), 
                       filters='ensembl_peptide_id', values=LRTdat2$Gene, mart=ensembl)
p0.05_go_info[p0.05_go_info$hgnc_symbol=='AGTRAP',]
p0.05_go_info[p0.05_go_info$hgnc_symbol=='C1QBP',]
p0.05_go_info[p0.05_go_info$hgnc_symbol=='PDIA3',]

p0.05_go_info[p0.05_go_info$namespace_1003=='biological_process',]

# check other datasets
ensembl_hs = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host='www.ensembl.org')  
ensembl_gg = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="ggallus_gene_ensembl", host='www.ensembl.org')  

p0.05_go_info_hs <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006', 'namespace_1003'), 
                       filters='hgnc_symbol', values=LRTdat2$hgnc_symbol, mart=ensembl_hs)
p0.05_go_info_gg <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006', 'namespace_1003'), 
                          filters='hgnc_symbol', values=LRTdat2$hgnc_symbol, mart=ensembl_gg)
p0.05_go_info_hs[p0.05_go_info_hs$hgnc_symbol=='AGTRAP',]
p0.05_go_info_hs[p0.05_go_info_hs$hgnc_symbol=='C1QBP',]
p0.05_go_info_hs[p0.05_go_info_hs$hgnc_symbol=='PDIA3',]
p0.05_go_info_hs[p0.05_go_info_hs$hgnc_symbol=='CCNI',]
p0.05_go_info_hs[p0.05_go_info_hs$hgnc_symbol=='MEA1',]
p0.05_go_info_hs[p0.05_go_info_hs$hgnc_symbol=='DNAJB2',]
p0.05_go_info_hs[p0.05_go_info_hs$hgnc_symbol=='MRPL42',]
p0.05_go_info_hs[p0.05_go_info_hs$hgnc_symbol=='GOT1',]
p0.05_go_info_hs[p0.05_go_info_hs$hgnc_symbol=='MRPS26',]
p0.05_go_info_hs[p0.05_go_info_hs$hgnc_symbol=='LYRM2',]

p0.05_go_info_gg[p0.05_go_info_gg$hgnc_symbol=='AGTRAP',]
p0.05_go_info_gg[p0.05_go_info_gg$hgnc_symbol=='C1QBP',]
p0.05_go_info_gg[p0.05_go_info_gg$hgnc_symbol=='PDIA3',]



# ----------------- search GO terms/info -----------------
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

alldata_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006', 'namespace_1003'), filters='ensembl_peptide_id', values=LRTdat$Gene, mart=ensembl)
head(alldata_go) #i don't think HBB is here because it doesn't have a zebra finch ensembl id

# let's try out a bar plot of name_1006, by LRT value, colored by namespace_1003
LRTdat3 <- merge(x=LRTdat[,c(1,6,8)], y=alldata_go, by.x='Gene', by.y='ensembl_peptide_id')
LRTdat3_MF <- LRTdat3[LRTdat3$namespace_1003 == 'molecular_function',]
ggplot(LRTdat3_MF, aes(x=LRT, y=name_1006, col=Significance_level)) + 
  geom_point() + theme_bw()

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

checkLRTdat <- LRTdat[LRTdat$signif_grp_001 == "signif",]

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

# ----------------- Non-metric multidimensional scaling -----------------
MF_NMDS = metaMDS(MFdat[13:ncol(MFdat)], k=2)
BP_NMDS = metaMDS(BPdat[13:ncol(BPdat)], k=2)
CC_NMDS = metaMDS(CCdat[13:ncol(CCdat)], k=2)
# not converging - doesn't matter if i add trymax=100, still doesn't converge 
# test different distance metrics

stressplot(MF_NMDS)
stressplot(BP_NMDS)
stressplot(CC_NMDS)

jpeg("all_NMDS.jpg", height=10, width=10, units="in", res=500)
par(mfrow=c(3,3))
ordiplot(MF_NMDS,type="n", main="MF, p=0.05")
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="red", lty=2)

ordiplot(MF_NMDS,type="n", main="MF, p=0.01")
ordihull(MF_NMDS, groups=MFdat$signif_grp_01, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(MF_NMDS, groups=MFdat$signif_grp_01, draw="polygon", label=T, show.groups="signif", col="red", lty=2)

ordiplot(MF_NMDS,type="n", main="MF, p=0.001")
ordihull(MF_NMDS, groups=MFdat$signif_grp_001, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(MF_NMDS, groups=MFdat$signif_grp_001, draw="polygon", label=T, show.groups="signif", col="red", lty=2)

ordiplot(BP_NMDS,type="n", main="BP, p=0.05")
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="red", lty=2)

ordiplot(BP_NMDS,type="n", main="BP, p=0.01")
ordihull(BP_NMDS, groups=BPdat$signif_grp_01, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(BP_NMDS, groups=BPdat$signif_grp_01, draw="polygon", label=T, show.groups="signif", col="red", lty=2)

ordiplot(BP_NMDS,type="n", main="BP, p=0.001")
ordihull(BP_NMDS, groups=BPdat$signif_grp_001, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(BP_NMDS, groups=BPdat$signif_grp_001, draw="polygon", label=T, show.groups="signif", col="red", lty=2)

ordiplot(CC_NMDS,type="n", main="CC, p=0.05")
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="red", lty=2)

ordiplot(CC_NMDS,type="n", main="CC, p=0.01")
ordihull(CC_NMDS, groups=CCdat$signif_grp_01, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(CC_NMDS, groups=CCdat$signif_grp_01, draw="polygon", label=T, show.groups="signif", col="red", lty=2)

ordiplot(CC_NMDS,type="n", main="CC, p=0.001")
ordihull(CC_NMDS, groups=CCdat$signif_grp_001, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(CC_NMDS, groups=CCdat$signif_grp_001, draw="polygon", label=T, show.groups="signif", col="red", lty=2)
dev.off()

jpeg("0.05_NMDS.jpg", height=4, width=6, units="in", res=500)
par(mfrow=c(1,3))
ordiplot(MF_NMDS,type="n", main="Molecular function")
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=F, show.groups="nonsig", col="darkgrey")
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=F, show.groups="signif", col="red", lty=2)

ordiplot(BP_NMDS,type="n", main="Biological process")
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=F, show.groups="nonsig", col="darkgrey")
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=F, show.groups="signif", col="red", lty=2)

ordiplot(CC_NMDS,type="n", main="Cellular component")
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=F, show.groups="nonsig", col="darkgrey")
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=F, show.groups="signif", col="red", lty=2)
dev.off()

par(mfrow=c(1,1))
# individual plots for p=0.05
jpeg("0.05MF_NMDS.jpg", height=10, width=10, units="in", res=500)
ordiplot(MF_NMDS,type="n", main="MF, p=0.05")
orditorp(MF_NMDS,display="species",cex=1.5,col="tomato",air=0.01)
orditorp(MF_NMDS,display="sites",cex=1,air=0.1)
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(MF_NMDS, groups=MFdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="red", lty=2)
#legend("topright", pch=c(16, 16), cex=c(1.5, 1.5), col=c("grey", "red"), c("Non-significant", "Significant"), bty="n")
dev.off()

jpeg("0.05BP_NMDS.jpg", height=10, width=10, units="in", res=500)
ordiplot(BP_NMDS,type="n", main="BP, p=0.05")
orditorp(BP_NMDS,display="species",cex=1.5,col="tomato",air=0.01)
orditorp(BP_NMDS,display="sites",cex=1,air=0.1)
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(BP_NMDS, groups=BPdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="red", lty=2)
#legend("topright", pch=c(16, 16), cex=c(1.5, 1.5), col=c("grey", "red"), c("Non-significant", "Significant"), bty="n")
dev.off()

jpeg("0.05CC_NMDS.jpg", height=10, width=10, units="in", res=500)
ordiplot(CC_NMDS,type="n", main="CC, p=0.05")
orditorp(CC_NMDS,display="species",cex=1.5,col="tomato",air=0.01)
orditorp(CC_NMDS,display="sites",cex=1,air=0.1)
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=T, show.groups="nonsig", col="darkgrey")
ordihull(CC_NMDS, groups=CCdat$signif_grp_05, draw="polygon", label=T, show.groups="signif", col="red", lty=2)
#legend("topright", pch=c(16, 16), cex=c(1.5, 1.5), col=c("grey", "red"), c("Non-significant", "Significant"), bty="n")
dev.off()

# ----------------- GO terms for each signif level -----------------
# new dataframes for each level of significance
p0.05 <- LRTdat[LRTdat$LRT >= 3.84,]
p0.01 <- LRTdat[LRTdat$LRT >= 6.64,]
p0.001 <- LRTdat[LRTdat$LRT >= 10.83,]
# read in list of genes not significant in codeml (at any level)
nonsig_genes <- LRTdat[LRTdat$LRT < 3.84,]

p0.05_go_names <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=p0.05$Gene, mart=ensembl)
dim(p0.05_go_names)
p0.05_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006', 'namespace_1003'), filters='ensembl_peptide_id', values=p0.05$Gene, mart=ensembl)
head(p0.05_go)
dim(p0.05_go)
p0.05_BP <- p0.05_go[p0.05_go$namespace_1003 == "biological_process",]
head(p0.05_BP)
dat_tab <- table(p0.05_BP$ensembl_peptide_id, p0.05_BP$name_1006) 
dat_tab[1:5, 1:5] 
#   #first column is an unidentified gene - that's why the colname is blank
#   #out of 80, there are only 60 with known biological_process functions
# write.csv(dat_tab, "dat_tab.csv")

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

# nonsig_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006', 'namespace_1003'), filters='ensembl_peptide_id', values=nonsig_genes$Gene, mart=ensembl)
# head(nonsig_go)
# dim(nonsig_go)
# nonsig_BP <- nonsig_go[nonsig_go$namespace_1003 == "biological_process",]
# nonsig_tab <- table(nonsig_BP$ensembl_peptide_id, nonsig_BP$name_1006)
# nonsig_tab[1:5, 1:5]










