# calculate LRT from codeml results, grab gene ontology categories
# null = Model A1 (NSsites = 2, model = 2, fix_omega = 1)
# test = Model A (NSsites = 2, model = 2, fix_omega = 0)

# For analysis where 1 high altitude species branch is tagged at a time
# 5 high alt species: Acast, Cviol, Mphoe, Pgiga, Ccoru

# 1. format lnL txt files and calculate LRT
# 2. calculate p-value and false discovery rate, extract gene name and function information
# 3. check codeml output omegas for genes with significant LRT - e.g., make sure the omega values make sense, filter out ones that are messed up
# 4. For remaining genes:
          # A. any overlaps between high alt species?
          # B. extract nucleotide site BEB scores

# Marisa Lim (c)2016

# load libraries
library(biomaRt)
library("GO.db")
library(qvalue)
library(ggplot2)

#library(vegan)
#library(dplyr)
#library(tidyr)
#library(ggrepel)

# set working directory
MLwd = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files/hi_alt_sp_lnL_analysis/"
setwd(MLwd)

# ----------------- 1. Calculate LRT -----------------

getLRT <- function(mysp, nuc_null, nuc_pos, mito_null, mito_pos){
  # nuclear dna files
  nu_nullmod <- read.table(nuc_null)
  nu_testmod <- read.table(nuc_pos)
  # give column names
  colnames(nu_nullmod) <- c("Gene", "lnL")
  colnames(nu_testmod) <- c("Gene", "lnL")
  # get gene name
  nu_replacebreak <- gsub(pattern="/", replacement=" ", x=nu_nullmod$Gene)
  nu_splitbreak <- strsplit(nu_replacebreak, split=" ")
  nu_genename1 <- sapply(nu_splitbreak, function(x){
    paste(x[[8]])
  })
  #head(nu_genename1)
  nu_replacebreak2 <- gsub(pattern="_", replacement=" ", x=nu_genename1)
  nu_splitbreak2 <- strsplit(nu_replacebreak2, split=" ")
  nu_genename <- sapply(nu_splitbreak2, function(x){
    paste(x[[1]])
  })
  #head(nu_genename)
  
  # mitochondrial dna files
  mt_nullmod <- read.table(mito_null)
  mt_testmod <- read.table(mito_pos)
  # give column names
  colnames(mt_nullmod) <- c("Gene", "lnL")
  colnames(mt_testmod) <- c("Gene", "lnL")
  # get gene name
  replacebreak <- gsub(pattern="/", replacement=" ", x=mt_nullmod$Gene)
  splitbreak <- strsplit(replacebreak, split=" ")
  genename1 <- sapply(splitbreak, function(x){
    paste(x[[8]])
  })
  #head(genename1)
  replacebreak2 <- gsub(pattern="_", replacement=" ", x=genename1)
  splitbreak2 <- strsplit(replacebreak2, split=" ")
  genename <- sapply(splitbreak2, function(x){
    paste(x[[1]])
  })
  #head(genename)
  
  # merge the information for both models 
  LRTdat_nu <- data.frame("Gene"=nu_genename, "Gene_nullmod"=nu_nullmod$Gene, "lnL_nullmod"=nu_nullmod$lnL, 
                          "Gene_testmod"=nu_testmod$Gene, "lnL_testmod"=nu_testmod$lnL)
  #head(LRTdat_nu)
  LRTdat_mito <- data.frame("Gene"=genename, "Gene_nullmod"=mt_nullmod$Gene, "lnL_nullmod"=mt_nullmod$lnL, 
                            "Gene_testmod"=mt_testmod$Gene, "lnL_testmod"=mt_testmod$lnL)
  #head(LRTdat_mito)
  
  # combine nuclear and mito datasets
  LRTdat <- rbind(LRTdat_nu, LRTdat_mito)

  # Want values > 0, ignore any < 0
  LRTdat$LRT <- 2*(LRTdat$lnL_testmod - LRTdat$lnL_nullmod)
  #head(LRTdat)
  
  LRTdat$Species <- mysp
  
  write.csv(LRTdat, paste(mysp, '_LRT.csv', sep=''))

}

getLRT('Acast', 'nu_Acast_lnL_null.txt', 'nu_Acast_lnL_pos.txt', 'Acast_lnL_null.txt', 'Acast_lnL_pos.txt')
getLRT('Cviol', 'nu_Cviol_lnL_null.txt', 'nu_Cviol_lnL_pos.txt', 'Cviol_lnL_null.txt', 'Cviol_lnL_pos.txt')
getLRT('Ccoru', 'nu_Ccoru_lnL_null.txt', 'nu_Ccoru_lnL_pos.txt', 'Ccoru_lnL_null.txt', 'Ccoru_lnL_pos.txt')
getLRT('Mphoe', 'nu_Mphoe_lnL_null.txt', 'nu_Mphoe_lnL_pos.txt', 'Mphoe_lnL_null.txt', 'Mphoe_lnL_pos.txt')
getLRT('Pgiga', 'nu_Pgiga_lnL_null.txt', 'nu_Pgiga_lnL_pos.txt', 'Pgiga_lnL_null.txt', 'Pgiga_lnL_pos.txt')

# ----------------- 2. Calculate p-values and FDR ---------------

mysp <- c('Acast', 'Cviol', 'Ccoru', 'Mphoe', 'Pgiga')
LRTdat <- c('Acast_LRT.csv', 'Cviol_LRT.csv', 'Ccoru_LRT.csv', 'Mphoe_LRT.csv', 'Pgiga_LRT.csv')

# Benjamini and Hochberg 1995 FDR method
getFDR_BH <- function(mysp, LRTdat){
  
  d <- data.frame('X'=1, 'Gene'='gene', 'Gene_nullmod'='genenullmod', 'lnL_nullmod'=1, 'Gene_testmod'='genetestmod', 'lnL_testmod'=1, 
                  'LRT'=1, 'Species'='species', 'p.value'=1, 'fdr_test'=1, 'FDR_ok'='yes', 'hgnc_symbol'='hgnc')
  
  for(i in 1:length(mysp)){
    
    theLRTdat <- read.csv(LRTdat[i], header=T)
    
    # extract p-values from LRT, use chi-square distribution
    thepvallist <- list()
    for(j in 1:nrow(theLRTdat)){
      theLRT <- theLRTdat$LRT[j]
      thepval <- pchisq(theLRT, df=1, lower.tail=FALSE)
      thepvallist[[j]] <- thepval
    }
    LRTdat_p <- data.frame(theLRTdat, "p-value"=unlist(thepvallist))
    #head(LRTdat_p)
    
    # FDR test Benjamini & Hochberg (1995) using full dataset 
    fdr_test <- p.adjust(p=LRTdat_p$p.value, method='BH')
    fdr_df <- as.data.frame(fdr_test)
    
    # this dataframe has uncorrected (p.value) and corrected (BH method, fdr_test) p-values
    # if you forget how to interpret, see this: https://www.r-bloggers.com/example-9-30-addressing-multiple-comparisons/
    LRTdat_fdr <- cbind(LRTdat_p, fdr_df) 
    
    # plot uncorrected and corrected p-values
    #ggplot(data=LRTdat_fdr) + 
    #  geom_line(aes(y=sort(p.value), x=X), col='blue', alpha=0.5, size=2) + 
    #  geom_line(aes(y=sort(fdr_test), x=X), col='tomato', alpha=0.5, size=2) +
    #  geom_hline(aes(yintercept=0.05), col='black', lty=2) +
    # geom_hline(aes(yintercept=0.01), col='grey', lty=2) +
    #  theme_bw() + ggtitle('P-values: uncorrected (blue) \nand BH corrected (red)')
    #ggsave(paste(mysp, '_pvals.jpg', sep=''), height=5, width=5, units='in', dpi=500)
    
    # subset data at significance level < 0.05
    LRTdat_fdr_ok <- LRTdat_fdr[LRTdat_fdr$fdr_test < 0.05,]
    LRTdat_fdr_ok$FDR_ok <- 'yes'
    
    # this is without the correction
    LRTdat_uncorrected <- LRTdat_fdr[LRTdat_fdr$p.value < 0.05,]
    LRTdat_uncorrected$FDR_ok <- 'no'
    
    # what's up with HBB? can add more candidate genes here later
    HBB <- LRTdat_fdr[LRTdat_fdr$Gene == 'HBB', ]
    HBB$FDR_ok <- 'na'
    
    # merge both dfs - I want to see which genes were significant before p-value correction
    LRTdatfinal <- rbind(LRTdat_uncorrected, LRTdat_fdr_ok, HBB)
    
    # now get gene names for the LRTdat_fdr_ok genes
    # if usual host not working today...25April16, the above host is an archived
    #ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="tguttata_gene_ensembl", host="dec2015.archive.ensembl.org")
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="tguttata_gene_ensembl", host='www.ensembl.org')  
    
    LRTdat_fdr_ok_withgenenames <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), 
                                         filters='ensembl_peptide_id', values=LRTdatfinal$Gene, mart=ensembl)
    LRTdat_checkomegas <- merge(x=LRTdatfinal, y=LRTdat_fdr_ok_withgenenames, by.x='Gene', by.y='ensembl_peptide_id', all.x=T)
    
    d <- rbind(d, LRTdat_checkomegas)
    
  }
  d2 <- d[-1,]
  
  write.csv(d2, 'gene_checkomegas_BHfdr.csv')
  
}
getFDR_BH(mysp, LRTdat)

# Storey and Tibshirani 2003 FDR method
getFDR_ST <- function(mysp, LRTdat){
  
  d <- data.frame('X'=1, 'Gene'='gene', 'Gene_nullmod'='genenullmod', 'lnL_nullmod'=1, 'Gene_testmod'='genetestmod', 'lnL_testmod'=1, 
                  'LRT'=1, 'Species'='species', 'p.value'=1, 'qvalues'=1, 'FDR_ok'='yes', 'hgnc_symbol'='hgnc')
  
  for(i in 1:length(mysp)){
    
    theLRTdat <- read.csv(LRTdat[i], header=T)
    
    # extract p-values from LRT, use chi-square distribution
    thepvallist <- list()
    for(j in 1:nrow(theLRTdat)){
      theLRT <- theLRTdat$LRT[j]
      thepval <- pchisq(theLRT, df=1, lower.tail=FALSE)
      thepvallist[[j]] <- thepval
    }
    LRTdat_p <- data.frame(theLRTdat, "p-value"=unlist(thepvallist))
    #head(LRTdat_p)
    
    # FDR test Storey and Tibshirani 2003 using full dataset 
    fdr_test <- qvalue(p=LRTdat_p$p.value)
    #summary(fdr_test)
    #hist(fdr_test)
    #plot(fdr_test)
    fdr_df <- data.frame(qvalues = fdr_test$qvalues)
    
    # this dataframe has uncorrected (p.value) and q-values (Storey and Tibshirani method, qvalues)
    # if you forget how to interpret, see this: https://www.r-bloggers.com/example-9-30-addressing-multiple-comparisons/
    LRTdat_fdr <- cbind(LRTdat_p, fdr_df) 
    
    # plot uncorrected and corrected p-values
    #ggplot(data=LRTdat_fdr) + 
    #  geom_line(aes(y=sort(p.value), x=X), col='blue', alpha=0.5, size=2) + 
    #  geom_line(aes(y=sort(qvalues), x=X), col='tomato', alpha=0.5, size=2) +
    #  geom_hline(aes(yintercept=0.05), col='black', lty=2) +
    # geom_hline(aes(yintercept=0.01), col='grey', lty=2) +
    #  theme_bw() + ggtitle('P-values: uncorrected (blue) \nand ST corrected (red)')
    #ggsave(paste(mysp, '_pvals.jpg', sep=''), height=5, width=5, units='in', dpi=500)
    
    # subset data at significance level < 0.05
    LRTdat_fdr_ok <- LRTdat_fdr[LRTdat_fdr$qvalues < 0.05,]
    LRTdat_fdr_ok$FDR_ok <- 'yes'
    
    # this is without the correction
    LRTdat_uncorrected <- LRTdat_fdr[LRTdat_fdr$p.value < 0.05,]
    LRTdat_uncorrected$FDR_ok <- 'no'
    
    # what's up with HBB? can add more candidate genes here later
    HBB <- LRTdat_fdr[LRTdat_fdr$Gene == 'HBB', ]
    HBB$FDR_ok <- 'na'
    
    # merge both dfs - I want to see which genes were significant before p-value correction
    LRTdatfinal <- rbind(LRTdat_uncorrected, LRTdat_fdr_ok, HBB)
    
    # now get gene names for the LRTdat_fdr_ok genes
    # if usual host not working today...25April16, the above host is an archived
    #ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="tguttata_gene_ensembl", host="dec2015.archive.ensembl.org")
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="tguttata_gene_ensembl", host='www.ensembl.org')  
    
    LRTdat_fdr_ok_withgenenames <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), 
                                         filters='ensembl_peptide_id', values=LRTdatfinal$Gene, mart=ensembl)
    LRTdat_checkomegas <- merge(x=LRTdatfinal, y=LRTdat_fdr_ok_withgenenames, by.x='Gene', by.y='ensembl_peptide_id', all.x=TRUE)
    
    d <- rbind(d, LRTdat_checkomegas)
    
  }
  d2 <- d[-1,]
  
  write.csv(d2, 'gene_checkomegas_STfdr.csv')
  
}
getFDR_ST(mysp, LRTdat)

# Note: the BH and ST methods for calculating correct p-value or qvalue yield same results.

# Look at gene_checkomegas.csv to check if codeml outputs for omega make sense

# Then, see if there are any overlaps in genes or functions between high alt species






# ----------------- search GO terms/info -----------------
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

alldata_genename <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), filters='ensembl_peptide_id', 
                          values=LRTdat$Gene, mart=ensembl)
head(alldata_genename)
dim(alldata_genename) # not including Hbb and HBD
write.csv(alldata_genename, 'Studygenes.csv')

alldata_go <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'name_1006', 'namespace_1003'), 
                    filters='ensembl_peptide_id', values=LRTdat$Gene, mart=ensembl)
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

# ----------------- Compare lit candidate vs. all my study genes list ------------
mygenelist <- read.csv('Studygenes.csv', header=T)
candlist <- read.csv('../Seq cap files/mycandgenes.csv', header=T)

head(mygenelist)
head(candlist)

length(mygenelist$hgnc_symbol)
length(candlist$Gene_name)

notinmystudy <- setdiff(candlist$Gene_name, mygenelist$hgnc_symbol)
length(notinmystudy)

inmystudy <- intersect(mygenelist$hgnc_symbol, candlist$Gene_name)
length(inmystudy)

length(candlist$Gene_name) - length(notinmystudy) 


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










