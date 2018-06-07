# Script for analyzing codeml results
# Marisa Lim (c)2018

# This version of the analysis was run on Stony Brook University's Seawulf cluster. 
# Batch and python scripts on Github: https://github.com/marisalim/Transcriptome_pipeline/tree/master/6Run_PAML

# ----- Project details ----------------

# For analyses comparing different branch tagging patterns:
# ---------------------------------------- #
# Tree      Foreground        Background
# ---------------------------------------- #
# B         High + mid        Low
# * (star)  High              Low + mid
# D         Low + mid         High
# E         Low               High + mid
# ---------------------------------------- #
# B vs. *: 2 branch-tagging patterns for testing for shared PSGs (positively selected genes) between high elevation species
  # varying where the mid-elevation species are grouped (either with high or low elev species)
# D vs. E: 2 branch-tagging patterns for testing for shared PSGs between low elevation species
  # varying where the mid-elevation species are grouped (either with high or low elev species)
# untested so far: tree A and C, which both exclude mid-elevation species from the tree topology entirely
  # not run yet because the analysis requires me to remove those species from the sequence alignments too
  # unless necessary, not doing yet, as i don't have an efficient method for removing species from these alignments easily
  # and not messing up the alignments for ~1000 genes without visually checking again in Geneious

# Objectives:
# 1. format lnL txt files and calculate LRT
  # null = Model A1 (NSsites = 2, model = 2, fix_omega = 1)
  # test = Model A (NSsites = 2, model = 2, fix_omega = 0)

# 2. calculate p-value and false discovery rate, extract gene name

# 3. check codeml output omegas for genes with significant LRT - 
  # e.g., make sure the omega values make sense, filter out ones that are messed up
  # Make a list of the genes you want to pull codeml output files for (to check omegas)
  
# 4. For remaining genes:
  # a. Are there overlapping functions?
  # b. Are there overlapping PSGs?
  # c. Are there overlapping sites?

# ----- Load libraries and setwd ----
library(ggplot2)
library(biomaRt) #biocLite('biomaRt')
library(ggrepel)

# add these later as needed (copied from previous script, so may not be needed here)
# library("GO.db") # biocLite("GO.db")
# library(qvalue) # biocLite('qvalue')
# library(biomartr)
# library(stringr)
# library(gridExtra)
# library(VennDiagram)
# library(ggtree)

# set working directory
MLwd = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files/Jan2018_codemlreanalysis/"
setwd(MLwd)

# ----- 1. Calculate LRT -----------------

getLRT <- function(mytree, nuc_null, nuc_pos, mito_null, mito_pos){
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
    paste(x[[2]]) #depending on file path structure, the indexing will change
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
    paste(x[[4]]) #depending on file path structure, the indexing will change
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
  #head(LRTdat_nu)

    LRTdat_mito <- data.frame("Gene"=genename, "Gene_nullmod"=mt_nullmod$Gene, "lnL_nullmod"=mt_nullmod$lnL, 
                            "Gene_testmod"=mt_testmod$Gene, "lnL_testmod"=mt_testmod$lnL)
  #head(LRTdat_mito)
  
  # combine nuclear and mito datasets
  LRTdat <- rbind(LRTdat_nu, LRTdat_mito)
  
  # Want values > 0, ignore any < 0
  LRTdat$LRT <- 2*(LRTdat$lnL_testmod - LRTdat$lnL_nullmod)
  #head(LRTdat)
  
  LRTdat$Tree <- mytree
  
  write.csv(LRTdat, paste(mytree, '_LRT.csv', sep=''))
  
}

getLRT('B', 'B_lnL_nuclnull.txt', 'B_lnL_nuclpos.txt', 'B_lnL_mitonull.txt', 'B_lnL_mitopos.txt')
getLRT('star', 'star_lnL_nuclnull.txt', 'star_lnL_nuclpos.txt', 'star_lnL_mitonull.txt', 'star_lnL_mitopos.txt')
getLRT('D', 'D_lnL_nuclnull.txt', 'D_lnL_nuclpos.txt', 'D_lnL_mitonull.txt', 'D_lnL_mitopos.txt')
getLRT('E', 'E_lnL_nuclnull.txt', 'E_lnL_nuclpos.txt', 'E_lnL_mitonull.txt', 'E_lnL_mitopos.txt')

# ----- 2. Calculate p-values and FDR ---------------

# Get p-values and do BH correction

# Benjamini and Hochberg 1995 FDR method
getFDR_BH <- function(mytree, myLRTdat){
  LRTdat <- read.csv(myLRTdat)
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
  ####################################################################################################################
  ### for this re-analysis, I am not sure this corrected p-value step is necessary, since I am
  ### currently considering each of these analyses as separate units, rather than combining their results... 
  ### 2/5/18: Liliana agrees that this isn't necessary
  ### I'll leave in the calculation however, in case this changes
  ####################################################################################################################
  fdr_test <- p.adjust(p=LRTdat_p$p.value, method='BH')
  fdr_df <- as.data.frame(fdr_test)
  
  # this dataframe has uncorrected (p.value) and corrected (BH method, fdr_test) p-values
  # if you forget how to interpret, see this: https://www.r-bloggers.com/example-9-30-addressing-multiple-comparisons/
  # and here: https://www.unc.edu/courses/2007spring/biol/145/001/docs/lectures/Nov12.html
  LRTdat_fdr <- cbind(LRTdat_p, fdr_df) 
  
  # plot uncorrected and corrected p-values
  ggplot(data=LRTdat_fdr) +
   geom_line(aes(y=sort(p.value), x=X), col='blue', alpha=0.5, size=2) +
   geom_line(aes(y=sort(fdr_test), x=X), col='tomato', alpha=0.5, size=2) +
   geom_hline(aes(yintercept=0.05), col='black', lty=2) +
   geom_hline(aes(yintercept=0.01), col='grey', lty=2) +
   theme_bw() + ylab('P-values: uncorrected (blue) \nand BH corrected (red)') + 
   ggtitle(paste(mytree, ' results'))
  ggsave(paste(mytree, '_pvals.jpg', sep=''), height=5, width=5, units='in', dpi=500)
  
  # subset data at significance level < 0.05
  LRTdat_fdr_ok <- LRTdat_fdr[LRTdat_fdr$fdr_test <= 0.05,]
  LRTdat_fdr_ok$FDR_ok <- 'yes, corrected p-value'
  
  # this is without the correction
  LRTdat_uncorrected <- LRTdat_fdr[LRTdat_fdr$p.value <= 0.05,]
  LRTdat_uncorrected$FDR_ok <- 'yes, uncorrected p-value'

  # merge both dfs - I want to see which genes were significant before p-value correction
  # this df ONLY has the PSGs that were detected either via p-value < 0.05 or corrected p-value < 0.05, 
  # genes that don't meet these significance levels are EXCLUDED
  LRTdatfinal <- rbind(LRTdat_uncorrected, LRTdat_fdr_ok)
  
  # now get gene names for the LRTdat_fdr_ok genes
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="tguttata_gene_ensembl", host='www.ensembl.org')  
  
  LRTdat_fdr_ok_withgenenames <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol'), 
                                       filters='ensembl_peptide_id', values=LRTdatfinal$Gene, mart=ensembl)
  LRTdat_checkomegas <- merge(x=LRTdatfinal, y=LRTdat_fdr_ok_withgenenames, by.x='Gene', by.y='ensembl_peptide_id', all.x=T)
  
  # need to fill in the hgnc column that i couldn't fill via biomart
  tguttata_annotations <- read.csv('../hi_alt_sp_lnL_analysis/tguttata_annots.csv')
  #head(tguttata_annotations)
  colnames(tguttata_annotations)[1] <- 'Gene'
  
  LRTdat_annotated <- merge(LRTdat_checkomegas, tguttata_annotations, by='Gene', all.x=TRUE)
  #head(LRTdat_annotated)
  #dim(LRTdat_annotated)
  
  # THIS IS THE FILE YOU WANT TO LOOK AT FOR PVALUE, FDR, ANNOTATIONS
  write.csv(LRTdat_annotated, paste('BHfdr_annot_TODOcheckomegas_', mytree, '.csv', sep=''))
  
  #plot the genes w/ and w/o correction (with Benjamini and Hochberg correction)
  #ggplot(LRTdat_annotated, aes(x=FDR_ok)) + geom_bar() + facet_wrap(~Tree, ncol=3) + theme_bw()
  ggplot(LRTdat_annotated, aes(x=fdr_test, y=p.value)) + geom_point() + facet_wrap(~FDR_ok, ncol=3) +
      geom_label_repel(aes(label=hgnc_symbol, fill=FDR_ok)) + theme_bw() +
      geom_vline(aes(xintercept=0.05), lty=2) +
      geom_hline(aes(yintercept=0.05), lty=2) +
      xlab('Benjamini Hochberg corrected p-value') + ylab('Uncorrected p-value') +
      ggtitle(paste(mytree, ' genes from paml analysis'))
    
  ggsave(paste(mytree, '_pvals-genes.jpg', sep=''), height=8, width=10, units='in', dpi=500)

}
getFDR_BH('B', 'B_LRT.csv')
getFDR_BH('star', 'star_LRT.csv')
getFDR_BH('D', 'D_LRT.csv')
getFDR_BH('E', 'E_LRT.csv')

# ----- 3. Make list of genes to pull output files for ----------------
# Make a text file list of the genes to check

genechecklist <- function(mytree){
  mydat <- read.csv(paste('BHfdr_annot_TODOcheckomegas_', mytree, '.csv', sep=''))
  # remove duplicate gene names (because it passed uncorrected and corrected p-value thresholds)
  mygenes <- unique(mydat$Gene)
  # add rest of path name to output file:
  myfullgenename <- list()
  for(i in 1:length(mygenes)){
    myfullgenename[[i]] <- paste('/gpfs/scratch/mclim/nucl_', mytree, '/nucl_outs', mytree, 
          '_all/', mygenes[i], '_', mytree, '_nucl_modApos.out', sep='')
  }
  
  write.table(unlist(myfullgenename), paste(mytree, '_codemlouts_tocheck.txt', sep=''), col.names=FALSE,
              row.names=FALSE, quote=FALSE)
}
genechecklist('B')
genechecklist('star')
genechecklist('D')
genechecklist('E')

# Next steps:
# 1. Just a few mito genes - manually edit the path
# 2. Have to edit the line endings because Windows makes it CRLF instead of LF
  # easy conversion: In Notepad++, go to search/replace window. Search for '\r\n' and replace with '\n' (no quotes in actual query)
# 3. Then use bash command to move those output files into a new dir 
  # (output files are on Seawulf (temp) and Drive B (permanent home))
  # use this: 
      # mkdir B_outputstocheck
      # xargs -a B_codemlouts_tocheck.txt cp -t ./B_outputstocheck/ 
# 4. Now the files to check (just the pos files right now) are in folders for each tree
  # check:
    # a. omega values
    # b. site class proportions (shouldn't be 0 if there is an omega value > 0...)
    # c. BEB sites
  # add results to 'Jan18_codemlout_checkBstar_DE.xlsx' 

# ----- 4. Analyze PSGs -----------
# 4a. Are there overlapping functions? 
  # For this set of analyses, look for overlapping functions between all the PSGs ID'd per tree
  # *** Pull GO terms, and any other pathway information
    # automating this with R packages hasn't worked very well since the output formats are not that useful or easy to interpret
    # probably best to do searches through the various gene/pathway databases manually


# 4b. Are there overlapping PSGs?
  # By definition, all PSGs per tree are shared by the tagged branch species
  # So, yes, there are overlapping PSGs
  
  # *** Nothing more to do here, except to find more annotation information for ensembl genes that I couldn't get annotated via biomart


# 4c. Are there overlapping sites? 
  # Positive sites fore foreground lineages were identified in several of the PSGs
  # These are recorded in 'Jan18_codemlout_checkBstar_DE.xlsx' column for each PSG

  # *** Next, if there is any information on the protein shape, check whether any of the sites are near interesting binding sites, etc. 
    # e.g., is there a biological reason why these particular sites might be evolving under positive selection







