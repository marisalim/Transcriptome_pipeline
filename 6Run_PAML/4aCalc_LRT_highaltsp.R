# calculate LRT from codeml results, grab gene ontology categories
# null = Model A1 (NSsites = 2, model = 2, fix_omega = 1)
# test = Model A (NSsites = 2, model = 2, fix_omega = 0)

# For analysis where 1 high altitude species branch is tagged at a time
# 5 high alt species: Acast, Cviol, Mphoe, Pgiga, Ccoru

# 1. format lnL txt files and calculate LRT
# 2. calculate p-value and false discovery rate, extract gene name
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
#library(VennDiagram)

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

# ----------------- 3. Look at gene_checkomegas.csv to check if codeml outputs for omega make sense ---------------
# There were several with proportion 2a and 2b = 0, discarded :(
# No overlap between species in genes.
# 9 genes, 2 for each species except only 1 for Pgigas




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












