# calculate LRT from codeml results, grab gene ontology categories
# null = Model A1 (NSsites = 2, model = 2, fix_omega = 1)
# test = Model A (NSsites = 2, model = 2, fix_omega = 0)

# For analysis where 1 species branch is tagged at a time for each of the 12 species

# 1. format lnL txt files and calculate LRT
# 2. calculate p-value and false discovery rate, extract gene name
# 3. check codeml output omegas for genes with significant LRT - e.g., make sure the omega values make sense, filter out ones that are messed up
# For remaining genes:
          # A. any overlaps between high alt species?
          # B. extract nucleotide site BEB scores
          # C. extract gene ontology (GO) and pathway information (see panther db)
# 4. Compare gene lists 
# 5. Plot GO and pathway information

# Marisa Lim (c)2016-17

# load libraries
library(biomaRt) #biocLite('biomaRt')
library("GO.db") # biocLite("GO.db")
library(qvalue) # biocLite('qvalue')
library(ggplot2)
library(ggrepel)
library(biomartr)
library(stringr)
library(gridExtra)
library(VennDiagram)
library(ggtree)

# set working directory
MLwd = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files/hi_alt_sp_lnL_analysis/"
setwd(MLwd)

# ----------------- 1. Calculate LRT -----------------
# P. malaris - 2255 wouldn't converge after many different starting values tried; remove this gene from all datasets
# there's also a redundant copy of HBB (as Hbb), remove Hbb and keep HBB

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
  dim(LRTdat_nu)
  # remove ENSTGUP00000002255 and Hbb
  LRTdat_nu2 <- LRTdat_nu[!LRTdat_nu$Gene == 'ENSTGUP00000002255', ]
  LRTdat_nu3 <- LRTdat_nu2[!LRTdat_nu2$Gene == 'Hbb', ]
  dim(LRTdat_nu2)
  dim(LRTdat_nu3)
  
  #head(LRTdat_nu)
  LRTdat_mito <- data.frame("Gene"=genename, "Gene_nullmod"=mt_nullmod$Gene, "lnL_nullmod"=mt_nullmod$lnL, 
                            "Gene_testmod"=mt_testmod$Gene, "lnL_testmod"=mt_testmod$lnL)
  #head(LRTdat_mito)
  
  # combine nuclear and mito datasets
  LRTdat <- rbind(LRTdat_nu3, LRTdat_mito)

  # Want values > 0, ignore any < 0
  LRTdat$LRT <- 2*(LRTdat$lnL_testmod - LRTdat$lnL_nullmod)
  #head(LRTdat)
  
  LRTdat$Species <- mysp
  
  write.csv(LRTdat, paste(mysp, '_LRT.csv', sep=''))

}

getLRT('Aamaz', 'nu_Aamaz_lnL_null.txt', 'nu_Aamaz_lnL_pos.txt', 'Aamaz_lnL_null.txt', 'Aamaz_lnL_pos.txt')
getLRT('Amela', 'nu_Amela_lnL_null.txt', 'nu_Amela_lnL_pos.txt', 'Amela_lnL_null.txt', 'Amela_lnL_pos.txt')
getLRT('Aviri', 'nu_Aviri_lnL_null.txt', 'nu_Aviri_lnL_pos.txt', 'Aviri_lnL_null.txt', 'Aviri_lnL_pos.txt')
getLRT('Ccoel', 'nu_Ccoel_lnL_null.txt', 'nu_Ccoel_lnL_pos.txt', 'Ccoel_lnL_null.txt', 'Ccoel_lnL_pos.txt')
getLRT('Cmuls', 'nu_Cmuls_lnL_null.txt', 'nu_Cmuls_lnL_pos.txt', 'Cmuls_lnL_null.txt', 'Cmuls_lnL_pos.txt')
getLRT('Phart', 'nu_Phart_lnL_null.txt', 'nu_Phart_lnL_pos.txt', 'Phart_lnL_null.txt', 'Phart_lnL_pos.txt')
getLRT('Pmala', 'nu_Pmala_lnL_null.txt', 'nu_Pmala_lnL_pos.txt', 'Pmala_lnL_null.txt', 'Pmala_lnL_pos.txt')
getLRT('Acast', 'nu_Acast_lnL_null.txt', 'nu_Acast_lnL_pos.txt', 'Acast_lnL_null.txt', 'Acast_lnL_pos.txt')
getLRT('Cviol', 'nu_Cviol_lnL_null.txt', 'nu_Cviol_lnL_pos.txt', 'Cviol_lnL_null.txt', 'Cviol_lnL_pos.txt')
getLRT('Ccoru', 'nu_Ccoru_lnL_null.txt', 'nu_Ccoru_lnL_pos.txt', 'Ccoru_lnL_null.txt', 'Ccoru_lnL_pos.txt')
getLRT('Mphoe', 'nu_Mphoe_lnL_null.txt', 'nu_Mphoe_lnL_pos.txt', 'Mphoe_lnL_null.txt', 'Mphoe_lnL_pos.txt')
getLRT('Pgiga', 'nu_Pgiga_lnL_null.txt', 'nu_Pgiga_lnL_pos.txt', 'Pgiga_lnL_null.txt', 'Pgiga_lnL_pos.txt')

# ----------------- 2. Calculate p-values and FDR ---------------

#1. need to read in all the LRT files and combine
multmerge = function(mypath){
  filenames=list.files(path=mypath, pattern='_LRT.csv', full.names=TRUE)
  datalist = lapply(filenames, function(x){
    read.csv(file=x, header=TRUE)
  })
  do.call(rbind, lapply(datalist, data.frame))
}
LRTdat_combined <- multmerge('./')
dim(LRTdat_combined)

# Benjamini and Hochberg 1995 FDR method
getFDR_BH <- function(LRTdat){
  
  # extract p-values from LRT, use chi-square distribution
  thepvallist <- list()
  for(j in 1:nrow(LRTdat)){
      theLRT <- LRTdat$LRT[j]
      thepval <- pchisq(theLRT, df=1, lower.tail=FALSE)
      thepvallist[[j]] <- thepval
    }
  LRTdat_p <- data.frame(LRTdat, "p-value"=unlist(thepvallist))
  #head(LRTdat_p)
    
    # FDR test Benjamini & Hochberg (1995) using full dataset 
  fdr_test <- p.adjust(p=LRTdat_p$p.value, method='BH')
  fdr_df <- as.data.frame(fdr_test)
    
    # this dataframe has uncorrected (p.value) and corrected (BH method, fdr_test) p-values
    # if you forget how to interpret, see this: https://www.r-bloggers.com/example-9-30-addressing-multiple-comparisons/
    # and here: https://www.unc.edu/courses/2007spring/biol/145/001/docs/lectures/Nov12.html
  LRTdat_fdr <- cbind(LRTdat_p, fdr_df) 
    
    # plot uncorrected and corrected p-values
    # ggplot(data=LRTdat_fdr) +
    #  geom_line(aes(y=sort(p.value), x=X), col='blue', alpha=0.5, size=2) +
    #  geom_line(aes(y=sort(fdr_test), x=X), col='tomato', alpha=0.5, size=2) +
    #  geom_hline(aes(yintercept=0.05), col='black', lty=2) +
    #  geom_hline(aes(yintercept=0.01), col='grey', lty=2) +
    #  theme_bw() + ylab('P-values: uncorrected (blue) \nand BH corrected (red)')
    # ggsave(paste(mysp, '_pvals.jpg', sep=''), height=5, width=5, units='in', dpi=500)
    
    # subset data at significance level < 0.05
  LRTdat_fdr_ok <- LRTdat_fdr[LRTdat_fdr$fdr_test <= 0.05,]
  LRTdat_fdr_ok$FDR_ok <- 'yes, corrected p-value'
    
    # this is without the correction
  LRTdat_uncorrected <- LRTdat_fdr[LRTdat_fdr$p.value <= 0.05,]
  LRTdat_uncorrected$FDR_ok <- 'yes, uncorrected p-value'
    
    # what's up with HBB? can add more candidate genes here later
    # HBB <- LRTdat_fdr[LRTdat_fdr$Gene == 'HBB', ]
    # HBB$FDR_ok <- 'check - Hbb'
    
    # merge both dfs - I want to see which genes were significant before p-value correction
  LRTdatfinal <- rbind(LRTdat_uncorrected, LRTdat_fdr_ok
                         # , HBB
                         )
    
    # now get gene names for the LRTdat_fdr_ok genes
    # if usual host not working today...25April16, the above host is an archived
    #ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="tguttata_gene_ensembl", host="dec2015.archive.ensembl.org")
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="tguttata_gene_ensembl", host='www.ensembl.org')  
    
  LRTdat_fdr_ok_withgenenames <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), 
                                         filters='ensembl_peptide_id', values=LRTdatfinal$Gene, mart=ensembl)
  LRTdat_checkomegas <- merge(x=LRTdatfinal, y=LRTdat_fdr_ok_withgenenames, by.x='Gene', by.y='ensembl_peptide_id', all.x=T)
    
  ggplot(LRTdat_checkomegas, aes(x=FDR_ok)) + geom_bar()

  write.csv(LRTdat_checkomegas, 'gene_checkomegas_BHfdr.csv')
  # write.csv(LRTdat_fdr, 'BHfdr_totalgenes.csv')
  
}
getFDR_BH(LRTdat_combined)

# need to fill in the hgnc column that i couldn't fill via biomart
tguttata_annotations <- read.csv('tguttata_annots.csv')
head(tguttata_annotations)
colnames(tguttata_annotations)[1] <- 'Gene'
BHdat <- read.csv('gene_checkomegas_BHfdr.csv')
head(BHdat)
dim(BHdat)
newdat <- merge(BHdat, tguttata_annotations, by='Gene', all.x=TRUE)
head(newdat)
dim(newdat)
#note: newdat has 1 extra entry compared to BHdat because ENSTGUP00000002606 has 2 annotations in tguttata
write.csv(newdat, 'checkomegas_BHfdr.csv')

#plot the genes w/ and w/o correction (with Benjamini and Hochberg correction)
ggplot(newdat, aes(x=FDR_ok)) + geom_bar() + facet_wrap(~Species, ncol=3) + theme_bw()
ggplot(newdat[newdat$FDR_ok == 'yes, corrected p-value',], aes(x=FDR_ok)) + geom_bar() +
  facet_wrap(~Species, ncol=4) + theme_bw()
for(i in 1:length(mysp)){
  thedat <- newdat[newdat$Species == mysp[i],]
  ggplot(thedat, aes(x=fdr_test, y=p.value)) + geom_point() + facet_wrap(~FDR_ok, ncol=3) +
    geom_label_repel(aes(label=hgnc_symbol, fill=FDR_ok)) + theme_bw() +
    geom_vline(aes(xintercept=0.05), lty=2) +
    geom_hline(aes(yintercept=0.05), lty=2) +
    xlab('Benjamini Hochberg corrected p-value') + ylab('Uncorrected p-value') +
    ggtitle(paste(mysp[i], ' genes from paml analysis'))

  ggsave(paste(mysp[i], '_pvals-genes.jpg', sep=''), height=8, width=10, units='in', dpi=500)

}

#plot corrected p-value in species pairs
#instead of plotting the dN/dS ratios, I can plot the p-values associated with 
#the significant or nonsignificant dN/dS for each gene
#obviously the point values are very different from dN/dS, but might be able to 
#interpret similarly
FULLcorrectedpvalsdat <- read.csv('BHfdr_totalgenes.csv')
head(FULLcorrectedpvalsdat)
unique(FULLcorrectedpvalsdat$Species)
ggplot(data=FULLcorrectedpvalsdat) + geom_bar(aes(x=fdr_test)) 

Aamaz_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Aamaz',]
Cmuls_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Cmuls',]
Acast_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Acast',]
Mphoe_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Mphoe',]
Amela_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Amela',]
Ccoru_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Ccoru',]
Phart_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Phart',]
Cviol_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Cviol',]
Pgiga_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Pgiga',]
Pmala_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Pmala',]
Ccoel_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Ccoel',]
Aviri_df <- FULLcorrectedpvalsdat[FULLcorrectedpvalsdat$Species=='Aviri',]

jpeg('FDRplots_speciespairs.jpg', height=7, width=7, units='in', res=600)
par(mfrow=c(3,3))
title('Plots of FDR for species pairs')
plot(-Cmuls_df$fdr_test, -Pgiga_df$fdr_test)
abline(b=1, a=0, lty=2)
abline(h=-0.05, col='red', lty=3)
abline(v=-0.05, col='red', lty=3)
plot(-Ccoel_df$fdr_test, -Cviol_df$fdr_test)
abline(b=1, a=0, lty=2)
abline(h=-0.05, col='red', lty=3)
abline(v=-0.05, col='red', lty=3)
plot(-Aviri_df$fdr_test, -Aamaz_df$fdr_test)
abline(b=1, a=0, lty=2)
abline(h=-0.05, col='red', lty=3)
abline(v=-0.05, col='red', lty=3)
plot(-Mphoe_df$fdr_test, -Amela_df$fdr_test)
abline(b=1, a=0, lty=2)
abline(h=-0.05, col='red', lty=3)
abline(v=-0.05, col='red', lty=3)
plot(-Phart_df$fdr_test, -Mphoe_df$fdr_test)
abline(b=1, a=0, lty=2)
abline(h=-0.05, col='red', lty=3)
abline(v=-0.05, col='red', lty=3)
plot(-Acast_df$fdr_test, -Ccoel_df$fdr_test)
abline(b=1, a=0, lty=2)
abline(h=-0.05, col='red', lty=3)
abline(v=-0.05, col='red', lty=3)
plot(-Pmala_df$fdr_test, -Ccoru_df$fdr_test)
abline(b=1, a=0, lty=2)
abline(h=-0.05, col='red', lty=3)
abline(v=-0.05, col='red', lty=3)
dev.off()

# **** Note: the BH and ST methods for calculating correct p-value or qvalue yield same results. ****
# Storey and Tibshirani 2003 FDR method
getFDR_ST <- function(mysp, LRTdat){

    # extract p-values from LRT, use chi-square distribution
    thepvallist <- list()
    for(j in 1:nrow(LRTdat)){
      theLRT <- LRTdat$LRT[j]
      thepval <- pchisq(theLRT, df=1, lower.tail=FALSE)
      thepvallist[[j]] <- thepval
    }
    LRTdat_p <- data.frame(LRTdat, "p-value"=unlist(thepvallist))
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

    write.csv(LRTdat_checkomegas, 'gene_checkomegas_STfdr.csv')

}
getFDR_ST(mysp, LRTdat)

# ----------------- 3. Look at gene_checkomegas.csv to check if codeml outputs for omega make sense ---------------
# There were several with proportion 2a and 2b = 0 >> discarded :(

# 2/17/17: I fixed FDR calculation, now need to compare old and new datasets for overlap in the codeml outputs I already checked
# correct file
correctdat <- read.csv('checkomegas_BHfdr.csv', header=TRUE)
head(correctdat)
dim(correctdat)
correctdat$matchby <- paste(correctdat$Species,'-', correctdat$Gene, '-', correctdat$FDR_ok, sep='')
# old file
omegacheckeddat <- read.csv('omega_check_9Feb17.csv')
head(omegacheckeddat)
names(omegacheckeddat)
dim(omegacheckeddat)
omegacheckeddat$matchby <- paste(omegacheckeddat$Species,'-', omegacheckeddat$Gene, '-', omegacheckeddat$FDR_ok, sep='')

# merge by codeml output path (has info about the gene and species)
test <- merge(x=omegacheckeddat, y=correctdat, by='matchby', all.y=TRUE); dim(test)
names(test)
# just keep relevant columns
test1 <- test[,c(1, 18:27, 13:15, 28:30)]; dim(test1)
# remove duplicate rows (these were introduced in the initial merge function I think)
test2 <- test1[!duplicated(test1),]; dim(test2)
# remove the genes that failed omega check
test3 <- test2[test2$Omega.check.short == 'Yes',]; dim(test3)

# there are 34 genes to look at functions, 145 without the FDR correction
# even before the p-value adjustment, there were a tiny proportion of genes 145/(12*947) * 100 = 1.2%
#   >> still lower than the assumed 5% false discovery rate??
# with correction 34/(12*947) * 100 = 0.3%
dim(test3[test3$FDR_ok.y == 'yes, corrected p-value',]) #34
dim(test3[test3$FDR_ok.y == 'yes, uncorrected p-value',]) #145
# this is a comparison of with and without correction for p-values
table(test2$FDR_ok.y, test2$Species.y, test2$Omega.check.short)
ggplot(test3, aes(x=FDR_ok.y)) + geom_bar() + facet_wrap(~Species.y, ncol=3) + theme_bw() 
# this is a comparison of genes under selection per species
test3a <- test3[test3$FDR_ok.y == 'yes, corrected p-value',]
ggplot(test3a, aes(x=Species.y)) + geom_bar(stat='count', fill='steelblue') + 
  theme_bw() + xlab('Species') + ylab('Number of genes under positive selection')
ggsave('NumgenesPosselbyspecies.jpg', height=5, width=5, units='in', dpi=600)

# does this have any phylogenetic pattern? look at on tree
# not really 



# plot the genes w/ and w/o correction (with Benjamini and Hochberg correction)
mysp <- c('Aamaz', 'Amela', 'Aviri', 'Ccoel', 'Cmuls', 'Phart', 'Pmala', 'Acast', 'Cviol', 'Ccoru', 'Mphoe', 'Pgiga')

for(i in 1:length(mysp)){
  thedat <- test3[test3$Species == mysp[i], ]
  ggplot(thedat, aes(x=fdr_test.y, y=p.value.y)) + geom_point() + 
    facet_wrap(~FDR_ok.y, ncol=3) + 
    geom_label_repel(aes(label=Associated.Gene.Name, fill=FDR_ok.y)) + theme_bw() + 
    geom_vline(aes(xintercept=0.05), lty=2) + 
    geom_hline(aes(yintercept=0.05), lty=2) +
    xlab('Benjamini Hochberg corrected p-value') + ylab('Uncorrected p-value') +
    ggtitle(paste(mysp[i], ' genes from paml analysis'))
  
  ggsave(paste(mysp[i], '_pvals-genesOMEGACHECKED.jpg', sep=''), height=8, width=10, units='in', dpi=500)
  
}

complete_fdr_dat <- read.csv('BHfdr_totalgenes.csv', header=T)
head(complete_fdr_dat)
notsigdat <- complete_fdr_dat[complete_fdr_dat$fdr_test > 0.05, ]
notsigdat$siglevel <- 'not significant'
sigdat <- complete_fdr_dat[complete_fdr_dat$fdr_test <= 0.05, ]
sigdat$siglevel <- 'significant'
complete_fdr_dat2 <- rbind(notsigdat, sigdat)
head(complete_fdr_dat2)

mycolors <- c('black', 'red')
ggplot(complete_fdr_dat2, aes(x=Gene, y=fdr_test, col=siglevel)) + 
  geom_point(alpha=0.5, size=3) + 
  scale_y_reverse() +
  scale_colour_manual(values=mycolors, guide=F) +
  # facet_wrap(~Species.y, ncol=3) + 
  geom_hline(aes(yintercept=0.05), col='tomato') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  ylab('Benjamini-Hochberg corrected p-values') +
  xlab('Genes')
ggsave('FDR_allgenes.jpg', height=6, width=6, units='in', dpi=600)
  
ggplot(complete_fdr_dat2, aes(x=fdr_test, fill=siglevel)) + 
  geom_histogram() + theme_bw()
  
write.csv(test3, 'addGOpath_tothisdataset_18feb17.csv')

# ----------------- 4. Compare lit candidate vs. all my study genes list ------------

# Genes identified under positive selection from my transcriptome study
my_pos_sel_list <- read.table('panther - all genes.txt')
# remove commas - 
my_pos_sel_list <- data.frame('gene'= gsub(',', '', my_pos_sel_list$V1))
# All genes analyzed in PAML in transcriptome study
transcriptome_genes_in_study <- read.csv('../Studygenes.csv', header=T)
colnames(transcriptome_genes_in_study) <- c('X', 'EnsemblProteinID', 'AssociatedGeneName')
# annotation is not so good here...
  # i manually added in the mitochondrial gene annotations (end of file), so that it will include COX1 and COX3



# Candidate genes from the literature (below, i redo this using the list from seq cap literature candidate genes)
litcandlist <- read.csv('../../Seq cap files/mycandgenes.csv', header=T)

head(my_pos_sel_list); dim(my_pos_sel_list)
head(transcriptome_genes_in_study); dim(transcriptome_genes_in_study)
head(litcandlist); dim(litcandlist)

# Q1: how many candidate genes from the lit were included in my transcriptome study?
notinmystudy <- setdiff(litcandlist$Gene_name, transcriptome_genes_in_study$AssociatedGeneName)
length(notinmystudy)

inmystudy <- intersect(transcriptome_genes_in_study$AssociatedGeneName, litcandlist$Gene_name)
length(inmystudy); length(litcandlist$Gene_name) - length(notinmystudy) 
length(inmystudy) - 5

# [1] "ISCA1"    "DNAJA1"   "CCNI"     "RPL37"    "AGTRAP"   "RAC3"     "DCXR"     "EPAS1"    "C1QBP"    "ATP6V1G3" "DDX5"     "MEA1"     "TUBA4A"  
# [14] "DNAJB2"   "RPL5"     "RHOA"     "HIF1AN"   "UBE2N"    "MRPL42"   "GOT1"     "RBX1"     "EGLN1"    "PSMD2"    "PSMC3"    "AP4S1"    "LYRM2"   
# [27] "CDKN1B"   "ACP1"     "COX1"     "COX3"  
# 30 previously identified candidate genes are included in my study
  # oops, a few of these are from my other paml analysis (when I tagged all high alt species together: 
    #AGTRAP, C1QBP, MYCBP, GOT1, ZNF638, MRPL42, MRPS26, LYRM2, PDIA3, DNAJB2, CCNI, MEA1)
    # but keep AGTRAP, DNAJB2, C1QBP ==> were identified by other studies!
    # SO, actually: 25 (take out CCNI, MEA1, LYRM2, GOT1, MRPL42)
        #doesn't include hemoglobin genes ==> BECAUSE THEY ARE NOT ANNOTATED IN input file
        # with HBB and HBAD: 25 + 2 = 27


 #Q2: how many of my positive selection genes overlap with previously id'd genes?
differentgenes <- setdiff(litcandlist$Gene_name, my_pos_sel_list$gene)
length(differentgenes)

samegenes <- intersect(litcandlist$Gene_name, my_pos_sel_list$gene)
length(samegenes) #DNAJA1

# a few differences might be due to inconsistent naming conventions, or gene names based on chicken rather than hgnc
# but overall, i'm seeing pretty different set of genes
# although their functions/pathways overlap

# 11/19/17: note this isn't complete either because it doesn't include the probe targets i made from transcriptome! :(
# let's try this again - match up candidate genes from lit (from seq capture pipeline, that were designed based on sequence in Anna's) 
# with ensembl IDs/annotations in zebra finch
# want to know what the ensembl IDs are so that I can see whether the genes were sequenced in my transcriptome dataset
#zebfinchensembl <- read.table('../Taeniopygia_guttata_gene_name.txt', header=T) # DOES NOT RUN: doesn't like the tabs or something?
#colnames(zebfinchensembl) <- c('EnsemblProteinID', 'AssociatedGeneName', 'WikiGeneDescription')
#zebfinchensembl <- zebfinchensembl[-1,]
#head(zebfinchensembl)
#dim(zebfinchensembl) #9724

# below is the list of candidate genes from lit that I found in Anna's hbird genome reference (used in seq capture pipeline)
# the original file was called /pylon5/bi4iflp/mlim/SeqCapData/mkref_Canna_targs/checkfastaheaders.txt
# 9/6/17: i think this does not include candidate genes I had from transcriptomes...
annacandsgenes <- read.table('../annacands_headers_seqcappipe.txt')
colnames(annacandsgenes) <- 'AssociatedGeneName'
head(annacandsgenes)
dim(annacandsgenes) #232
annacandsgenes$AssociatedGeneName <- as.character(annacandsgenes$AssociatedGeneName)
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q68Y81_CHICK'] <- 'FANCD2'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'C16orf75'] <- 'RMI2'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'NOS2_CHICK'] <- 'NOS2'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q5F4B9_CHICK'] <- 'ABCD3'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'HDAC9_CHICK'] <- 'HDAC9'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'LAMB1_CHICK'] <- 'LAMB1'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q58I02_CHICK'] <- 'PBEF1'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q788Q2_CHICK'] <- 'HGF'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q5ZMU9_CHICK'] <- 'VCP'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q5ZL56_CHICK'] <- 'ACADS'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'C15orf42'] <- 'TICRR'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'ACE_CHICK'] <- 'ACE'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'UTP15_CHICK'] <- 'UTP15'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q5QHR9_CHICK'] <- 'ADAM17'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'NRG1_CHICK'] <- 'NRG1'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'KS6A5_CHICK'] <- 'RPS6KA5'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'NTRK2_CHICK'] <- 'NTRK2'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'RS6_CHICK'] <- 'RPS6'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q98SN3_CHICK'] <- 'ARNT'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'PDIA3_CHICK'] <- 'PDIA3'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'CO6A1_CHICK'] <- 'COL6A1'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'CAD13_CHICK'] <- 'CDH13'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'THB_CHICK'] <- 'THRB'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q9DGL8_CHICK'] <- 'JSC'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'C17orf60'] <- 'MILR1'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q2ACD4_CHICK'] <- 'MRPS26'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'Q5ZJW2_CHICK'] <- 'NDUFB8'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'VINC_CHICK'] <- 'VCL'
annacandsgenes$AssociatedGeneName[annacandsgenes$AssociatedGeneName == 'KCMA1_CHICK'] <- 'KCNMA1'
annacandsgenes$AssociatedGeneName <- as.factor(annacandsgenes$AssociatedGeneName)

# how many cand genes are in my study already?
inzeb <- merge(annacandsgenes, transcriptome_genes_in_study, by='AssociatedGeneName')
head(inzeb)
dim(inzeb) #28 (add Hbb to this, so 29 genes included). BUT, seems to be some missing ...e.g., EGLN1
# *** missing genes (e.g., EGLN1) because they were designed from transcripts not anna's genome. ***

inzeb2 <- inzeb
inzeb2$AssociatedGeneName <- as.character(inzeb2$AssociatedGeneName)
inzeb2$AssociatedGeneName[inzeb2$AssociatedGeneName == 'ACADS'] <- 'Q5ZL56_CHICK'
inzeb2$AssociatedGeneName[inzeb2$AssociatedGeneName == 'MRPS26'] <- 'Q2ACD4_CHICK'
inzeb2$AssociatedGeneName[inzeb2$AssociatedGeneName == 'RPS6'] <- 'RS6_CHICK'
inzeb2$AssociatedGeneName[inzeb2$AssociatedGeneName == 'PDIA3'] <- 'PDIA3_CHICK'
inzeb2$AssociatedGeneName <- as.factor(inzeb2$AssociatedGeneName)

# candidate genes to blast search for against Trinity transcriptome assemblies:
toblast <- setdiff(annacandsgenes$AssociatedGeneName, inzeb$AssociatedGeneName)
length(toblast) #203
toblastdf <- data.frame('Gene'=toblast)
#have to change gene name back (so that it will match name in ref fasta file)
toblastdf$Gene <- as.character(toblastdf$Gene)
toblastdf$Gene[toblastdf$Gene == 'FANCD2'] <- 'Q68Y81_CHICK'
toblastdf$Gene[toblastdf$Gene == 'RMI2'] <- 'C16orf75'
toblastdf$Gene[toblastdf$Gene == 'NOS2'] <- 'NOS2_CHICK'
toblastdf$Gene[toblastdf$Gene == 'ABCD3'] <- 'Q5F4B9_CHICK'
toblastdf$Gene[toblastdf$Gene == 'HDAC9'] <- 'HDAC9_CHICK'
toblastdf$Gene[toblastdf$Gene == 'LAMB1'] <- 'LAMB1_CHICK'
toblastdf$Gene[toblastdf$Gene == 'PBEF1'] <- 'Q58I02_CHICK'
toblastdf$Gene[toblastdf$Gene == 'HGF'] <- 'Q788Q2_CHICK'
toblastdf$Gene[toblastdf$Gene == 'VCP'] <- 'Q5ZMU9_CHICK'
toblastdf$Gene[toblastdf$Gene == 'ACADS'] <- 'Q5ZL56_CHICK'
toblastdf$Gene[toblastdf$Gene == 'TICRR'] <- 'C15orf42'
toblastdf$Gene[toblastdf$Gene == 'ACE'] <- 'ACE_CHICK'
toblastdf$Gene[toblastdf$Gene == 'UTP15'] <- 'UTP15_CHICK'
toblastdf$Gene[toblastdf$Gene == 'ADAM17'] <- 'Q5QHR9_CHICK'
toblastdf$Gene[toblastdf$Gene == 'NRG1'] <- 'NRG1_CHICK'
toblastdf$Gene[toblastdf$Gene == 'RPS6KA5'] <- 'KS6A5_CHICK'
toblastdf$Gene[toblastdf$Gene == 'NTRK2'] <- 'NTRK2_CHICK'
toblastdf$Gene[toblastdf$Gene == 'RPS6'] <- 'RS6_CHICK'
toblastdf$Gene[toblastdf$Gene == 'ARNT'] <- 'Q98SN3_CHICK'
toblastdf$Gene[toblastdf$Gene == 'PDIA3'] <- 'PDIA3_CHICK'
toblastdf$Gene[toblastdf$Gene == 'COL6A1'] <- 'CO6A1_CHICK'
toblastdf$Gene[toblastdf$Gene == 'CDH13'] <- 'CAD13_CHICK'
toblastdf$Gene[toblastdf$Gene == 'THRB'] <- 'THB_CHICK'
toblastdf$Gene[toblastdf$Gene == 'JSC'] <- 'Q9DGL8_CHICK'
toblastdf$Gene[toblastdf$Gene == 'MILR1'] <- 'C17orf60'
toblastdf$Gene[toblastdf$Gene == 'MRPS26'] <- 'Q2ACD4_CHICK'
toblastdf$Gene[toblastdf$Gene == 'NDUFB8'] <- 'Q5ZJW2_CHICK'
toblastdf$Gene[toblastdf$Gene == 'VCL'] <- 'VINC_CHICK'
toblastdf$Gene[toblastdf$Gene == 'KCNMA1'] <- 'KCMA1_CHICK'
toblastdf$Gene <- as.factor(toblastdf$Gene)

write.table(toblastdf, '15may17_candgenestoblast_toTrinityassemblies.txt', row.names=FALSE)
write.table(inzeb2, '15may17_candgenestoEXCLUDEfromblast.txt', row.names=FALSE) #these are the sequences to remove from fasta ref file for blast search
# 9/6/17 note: see parseblastresults.r for rest of this part - wasn't able to get overlaps for all 12 species

# ----------------- 5. Plot GO and pathway information ---------------

# ------ 5a. Exploring the GO data ------
# Genes identified under positive selection from my transcriptome study
my_pos_sel_list <- read.table('panther - all genes.txt')
# remove commas - 
my_pos_sel_list2 <- list(my_pos_sel_list$V1)
my_pos_sel_list <- gsub(',', '', my_pos_sel_list$V1)
my_pos_sel_df <- data.frame(my_pos_sel_list)
colnames(my_pos_sel_df) <- 'gene'

mygoslim <- getGO(organism = 'Homo sapiens', genes = c('CACTIN','TMEM38A','RRP8','ANXA6','TCEA3','DNAJA1',
                                                       'YRDC','MDH1','GPD1','PSAP',
                                                       'TPM1','ATG9A','CSTL1','NDUFB10','NOB1','MFGE8',  
                                                       'SDHA','TIMM21','GPI','EIF4A2','MLF1','ATP6V1D','JPH1','UBXN4',
                                                       'THRSP','NDUFB4','CLPB','FEM1A','ND1','COX1'), 
                  filters= 'hgnc_symbol')
mygoslim
dim(mygoslim)
### the go slim categories are not in the same order, but everything that follows MF, BP, or CC will be in that category

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host='www.ensembl.org')  
ensembl@attributes$name
GOpath_info <- getBM(attributes=c('hgnc_symbol','name_1006', 'namespace_1003'), 
                     filters='hgnc_symbol', values=my_pos_sel_df$gene, mart=ensembl)
# add attribute 'definition_1006' for more details
dim(GOpath_info)
GOpath_info[GOpath_info$name_1006 == '',] <- 'unknown'
ggplot(GOpath_info[GOpath_info$namespace_1003 == 'biological_process',], aes(x=name_1006)) + 
  geom_bar() + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggplot(GOpath_info[GOpath_info$namespace_1003 == 'biological_process',], aes(x=goslim_goa_description)) +
  geom_bar() +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))

GOpath_info_repeatedonly <- GOpath_info[duplicated(GOpath_info$name_1006),]; dim(GOpath_info_repeatedonly)
ggplot(GOpath_info_repeatedonly[GOpath_info_repeatedonly$namespace_1003 == 'biological_process',], 
       aes(x=name_1006)) + geom_bar() + theme_bw() +
  theme(axis.text.x=element_text(angle=55, hjust=1))

# ------ 5b. Plotting the GO data ------
# read in GO data (manually formatted) from Panther -- these are the complete GO lists, not GO slim (which were incomplete from panther)
mygodat <- read.csv('GOterms_forwordcloud_bar.csv', head=T)
mygodat_mf <- mygodat[mygodat$GO.category == 'MF',]
mygodat_bp <- mygodat[mygodat$GO.category == 'BP',]
mygodat_cc <- mygodat[mygodat$GO.category == 'CC',]

# makes new column that wraps text so it's not so long
mygodat_bp$GOdescription.wrap <- str_wrap(mygodat_bp$GO.database.description, width=30)

mycolors <- c("#c51b7d", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#a6cee3", "#fdbf6f", "#ff7f00", "#cab2d6")
ggplot(mygodat_bp[mygodat_bp$Elev.2000m == 'high',], aes(x=GO.database.description)) +
  geom_bar(position='stack', aes(fill=Species)) + coord_flip() +
  # facet_wrap(~Species) + 
  theme_bw() + scale_fill_manual(values=mycolors) + ylab('Count') +
  # theme(axis.text.x=element_text(angle=55, hjust=1)) + 
  xlab('Biological Process Descriptions')
ggsave('BP_byhigh2000mspecies.jpg', width=10, height=10, units='in', dpi=600)





# this is still quite messy to look at
# how about a venn diagram between elevation categories?
# ------ 5b.(i) Biological process subsets -----
nrow(subset(mygodat_bp, Elev.2000m == 'high')) #108
nrow(subset(mygodat_bp, Elev.2000m == 'low')) #116
nrow(subset(mygodat_bp, Elev.hbb == 'high')) #66
nrow(subset(mygodat_bp, Elev.hbb == 'low')) #158

overlapElev.2000 <- intersect(mygodat_bp[mygodat_bp$Elev.2000m == 'high',]$GO.database.description, 
                              mygodat_bp[mygodat_bp$Elev.2000m == 'low', ]$GO.database.description)
length(overlapElev.2000) #9

overlapElev.hbb <- intersect(mygodat_bp[mygodat_bp$Elev.hbb == 'high',]$GO.database.description, 
                              mygodat_bp[mygodat_bp$Elev.hbb == 'low', ]$GO.database.description)
length(overlapElev.hbb) #9

# plot venn diagram
grid.newpage() # Elev.2000m
bpplot.2000 <- draw.pairwise.venn(area1 = 108, area2 = 116, cross.area = 9, 
                   category = c("High altitude species", "Low altitude species"),
                   fill=c('blue', 'tomato'), lty=rep('blank',2),
                   cat.pos = c(0,0))

grid.newpage() #Elev.hbb
bpplot.hbb <- draw.pairwise.venn(area1 = 66, area2 = 91, cross.area = 9, 
                   category = c("High altitude species", "Low altitude species"),
                   fill=c('blue', 'tomato'), lty=rep('blank',2),
                   cat.pos = c(0,0))

jpeg(filename = "GOtermvenn_bp2000m.jpeg", height=5, width=5, units='in', res=500)
grid.draw(bpplot.2000)
dev.off()
jpeg(filename = "GOtermvenn_bphbb.jpeg", height=5, width=5, units='in', res=500)
grid.draw(bpplot.hbb)
dev.off()
# ------ 5b.(ii) Molecular function subsets -----
nrow(subset(mygodat_mf, Elev.2000m == 'high')) #55
nrow(subset(mygodat_mf, Elev.2000m == 'low')) #66
nrow(subset(mygodat_mf, Elev.hbb == 'high')) #30
nrow(subset(mygodat_mf, Elev.hbb == 'low')) #91

overlapElev.2000 <- intersect(mygodat_mf[mygodat_mf$Elev.2000m == 'high',]$GO.database.description, 
                              mygodat_mf[mygodat_mf$Elev.2000m == 'low', ]$GO.database.description)
length(overlapElev.2000) #11

overlapElev.hbb <- intersect(mygodat_mf[mygodat_mf$Elev.hbb == 'high',]$GO.database.description, 
                             mygodat_mf[mygodat_mf$Elev.hbb == 'low', ]$GO.database.description)
length(overlapElev.hbb) #12

# plot venn diagram
grid.newpage() # Elev.2000m
mfplot.2000 <- draw.pairwise.venn(area1 = 55, area2 = 66, cross.area = 11, 
                   category = c("High altitude species", "Low altitude species"),
                   fill=c('blue', 'tomato'), lty=rep('blank',2),
                   cat.pos = c(0,0))

grid.newpage() #Elev.hbb
mfplot.hbb <- draw.pairwise.venn(area1 = 30, area2 = 91, cross.area = 12, 
                   category = c("High altitude species", "Low altitude species"),
                   fill=c('blue', 'tomato'), lty=rep('blank',2),
                   cat.pos = c(0,0))

jpeg(filename = "GOtermvenn_mf2000m.jpeg", height=5, width=5, units='in', res=500)
grid.draw(mfplot.2000)
dev.off()
jpeg(filename = "GOtermvenn_mfhbb.jpeg", height=5, width=5, units='in', res=500)
grid.draw(mfplot.hbb)
dev.off()
# ------ 5b.(iii) Cellular component subsets -----
nrow(subset(mygodat_cc, Elev.2000m == 'high')) #96
nrow(subset(mygodat_cc, Elev.2000m == 'low')) #77
nrow(subset(mygodat_cc, Elev.hbb == 'high')) #68
nrow(subset(mygodat_cc, Elev.hbb == 'low')) #105

overlapElev.2000 <- intersect(mygodat_cc[mygodat_cc$Elev.2000m == 'high',]$GO.database.description, 
                              mygodat_cc[mygodat_cc$Elev.2000m == 'low', ]$GO.database.description)
length(overlapElev.2000) #15

overlapElev.hbb <- intersect(mygodat_cc[mygodat_cc$Elev.hbb == 'high',]$GO.database.description, 
                             mygodat_cc[mygodat_cc$Elev.hbb == 'low', ]$GO.database.description)
length(overlapElev.hbb) #6

# plot venn diagram
grid.newpage() # Elev.2000m
ccplot.2000 <- draw.pairwise.venn(area1 = 96, area2 = 77, cross.area = 15, 
                   category = c("High altitude species", "Low altitude species"),
                   fill=c('blue', 'tomato'), lty=rep('blank',2),
                   cat.pos = c(0,0))

grid.newpage() #Elev.hbb
ccplot.hbb <- draw.pairwise.venn(area1 = 68, area2 = 105, cross.area = 16, 
                   category = c("High altitude species", "Low altitude species"),
                   fill=c('blue', 'tomato'), lty=rep('blank',2),
                   cat.pos = c(0,0))

jpeg(filename = "GOtermvenn_cc2000m.jpeg", height=5, width=5, units='in', res=500)
grid.draw(ccplot.2000)
dev.off()
jpeg(filename = "GOtermvenn_cchbb.jpeg", height=5, width=5, units='in', res=500)
grid.draw(ccplot.hbb)
dev.off()






# ------ 5c. Plotting BP GO terms from Panther low vs. high
bpgo_dat <- read.csv('../transcriptome paper drafts/PantherBPGO_plotexplore.csv')
head(bpgo_dat)

ggplot(data=bpgo_dat) + geom_bar(aes(x=GO.database.BP.complete)) +
  facet_wrap(~Elevation) + theme_bw() + coord_flip() + xlab('Biological Process GO term')

ggplot(data=bpgo_dat) + geom_bar(aes(x=GO.database.BP.complete, fill=Elevation), position='stack') +
  theme_bw() + coord_flip() + xlab('Biological Process GO term')


