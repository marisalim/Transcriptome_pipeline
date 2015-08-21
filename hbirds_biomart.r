# Marisa Lim (c)2015
MLwd = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files"
setwd(MLwd)
# -------------------------------------
# find RBHs shared between all species
# -------------------------------------
rbh <- read.csv("rbh_all.csv", header=T)
head(rbh)
dim(rbh)

vecA <- as.vector(rbh$A)
vecB <- as.vector(rbh$B)
vecC <- as.vector(rbh$C)
vecD <- as.vector(rbh$D)
vecE <- as.vector(rbh$E)
vecF <- as.vector(rbh$F)
vecG <- as.vector(rbh$G)
vecH <- as.vector(rbh$H)
vecI <- as.vector(rbh$I)
vecJ <- as.vector(rbh$J)
vecK <- as.vector(rbh$K)
vecL <- as.vector(rbh$L)
hits <- Reduce(intersect, list(vecA, vecB, vecC, vecD, vecE, vecF, vecG, vecH, vecI, vecJ, vecK, vecL))
length(hits)

head(hits)
hitsvec <- as.vector(hits)
length(hitsvec) 

# ---------------------------------------------------------
# Use biomaRt to find gene names and GO info for these hits
# ---------------------------------------------------------
library(biomaRt)
#biocLite("GO.db")
library("GO.db")

ensembl = useMart("ensembl", dataset="tguttata_gene_ensembl")
filters = listFilters(ensembl)
head(filters)
attributes = listAttributes(ensembl)
head(attributes)

# test code for biomaRt
grep("description",attributes$name, ignore.case=T, value=T)

test <- 'ENSTGUP00000000081'
getBM(attributes=c('ensembl_peptide_id', 'go_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=test, mart=ensembl)

go_search <- getBM(attributes="go_id", filters="ensembl_peptide_id", values = test, mart = ensembl)
Term(go_search$go_id)

as.data.frame(Term(go_search$go_id))

# look at descriptions for the first 5 go terms
for(go in go_search$go_id[1:5]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

#test with vector of peptides
test2 <- c('ENSTGUP00000012872', 'ENSTGUP00000000081', 'ENSTGUP00000010502', 'ENSTGUP00000006586', 'ENSTGUP00000003270',
           'ENSTGUP00000005098', 'ENSTGUP00000014035', 'ENSTGUP00000006602', 'ENSTGUP00000001084', 'ENSTGUP00000001948')
test2go <- getBM(attributes=c('ensembl_peptide_id', 'go_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=test2, mart=ensembl)
Term(test2go$go_id)
as.data.frame(Term(test2go$go_id))
#one of the GO terms is obsolete, use na.omit to circumvent
na.omit(as.data.frame(Term(test2go$go_id)))

# now with real data:
hits_go <- getBM(attributes=c('ensembl_peptide_id', 'go_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=hitsvec, mart=ensembl)
head(hits_go)
dim(hits_go)
class(hits_go)

goterms <- as.data.frame(Term(hits_go$go_id))
head(goterms)
class(goterms)

hits_go2 <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'description'), filters='ensembl_peptide_id', values=hitsvec, mart=ensembl)
head(hits_go2)
dim(hits_go2)

write.csv(hits_go2, "Sharedhits_gene_descript.csv")

hits_list <- data.frame(Ensembl_id = hits_go2[,1])
write.csv(hits_list, "Sharedhits_list.csv")

# -------------------------------------------------------------------------
# Match ensembl id of shared hits to contig ID (need this to find sequence)
# Get sequences from matched contig ids
# -------------------------------------------------------------------------
require(plyr)
library(seqinr)

# Function to match shared hits list to contig ids, and then use contig ids to pull out the sequences [note: the rbhhits and rhbseqs files need to be made first]
matchcontig_getseq <- function(){
  # load shared genes list
  hits_list2 <- read.csv("Sharedhits_list.csv", header=T)
  hits_list2 <- data.frame("Ensembl_id"=hits_list2$Ensembl_id)
  #head(hits_list2)
  
  # Create infile lists
  rbhfolder = "RBH/"
  rbhfiles <- list.files(paste(rbhfolder,sep="/"),pattern=".out",full.name=T,recursive=T)
  fastafolder = "trinity_assembly/"
  fastafiles <- list.files(paste(fastafolder,sep="/"),pattern=".fasta",full.name=T,recursive=T)

  # Create sample id vector
  sampleids <- c("a_Phart", "b_Cmuls", "c_Pmala", "d_Acast", "e_Ccoru", "f_Ccoel", "g_Cviol", "h_Aamaz", "i_Mphoe", "j_Aviri", "k_Amela", "l_Pgiga")
  
  # Match shared hits to species rbh
  for(i in 1:length(rbhfiles)){
    # load infile_rbh
    infile_rbh <- read.table(rbhfiles[i])
    colnames(infile_rbh) <- c("contig_id", "Ensembl_id")
    
    # get sample id
    samp <- sampleids[i]
    
    # add samp column to infile_rbh
    infile_rbh["Sample_id"] <- NA
    infile_rbh$Sample_id <- samp
    # match contig id from species rbh to shared gene list
    hits2 <- join(hits_list2, infile_rbh, by="Ensembl_id")
    # modify contig id so that it matches name in assembly fasta file
    replacebreak <- gsub(pattern="\\|", replacement=" ", x=hits2$contig_id)
    splitbreak <- strsplit(replacebreak, split=" ")
    contigname <- sapply(splitbreak, function(x){
      paste(x[[1]])
    })
    # add modified contig id to df
    hits3 <- data.frame("Ensembl_id"=hits2$Ensembl_id, "Contig_id_long"=hits2$contig_id, "Contig_id"=contigname, "Sample_id"=hits2$Sample_id)
    write.table(hits3, paste("rbhhits/",samp,"_hits3.txt", sep=""))
  }
  
  # Need to load hits3 files for next loop
  hits3folder <- "rbhhits/"
  hits3files <- list.files(paste(hits3folder, sep="/"), pattern=".txt", full.name=T, recursive=T)
  
  # Match contig ID to get sequences from fasta files
  for(j in 1:length(fastafiles)){
    # load infile_fasta
    infile_fasta <- read.fasta(fastafiles[j], seqtype="DNA", as.string=T, forceDNAtolower=F)
    
    # make fasta list to df
    fastadf <- ldply(.data=infile_fasta)
    colnames(fastadf) <- c("Contig_id", "Sequence")
    
    # load the hits3 lists
    hits3 <- read.table(hits3files[j])
    # get sample id
    samp <- sampleids[j]
    # match contig id from rbh shared hits for given species to contig id in assembly file
    seqs <- join(hits3, fastadf, by="Contig_id") 
    # save file
    write.table(seqs, paste("rbhseqs/",samp,"_seqs.txt", sep="")) 
  }
}
matchcontig_getseq()

# ---------------------------------------------------------
# Match up the sequences for each gene from the 12 species
# ---------------------------------------------------------
require(plyr)
require(seqinr)

# load shared genes list
hits_list2 <- read.csv("Sharedhits_list.csv", header=T)
hits_list2 <- data.frame("Ensembl_id"=hits_list2$Ensembl_id)
head(hits_list2)

# TODO: can probably wrap this into the function
# for each of the 12 seq files, grab all rows that match given ensembl id from hits_list2
seqfolder = "rbhseqs/"
seqfiles <- list.files(paste(seqfolder,sep="/"),pattern=".txt",full.name=T,recursive=T)
myseq <- list()
for(i in 1:nrow(hits_list2)){
  # define ensembl id to find
  target <- hits_list2[i,]
  # find target in seq files
  for(j in 1:length(seqfiles)){
    # define seq file
    seqfile <- read.table(seqfiles[j])
    myseq[[j]] <- seqfile[seqfile$Ensembl_id == target, ]
  }
  
}



for(i in 1:2){
  # define ensembl id to find
  target <- hits_list2[i,]
  # find target in seq files, store in myseq list
  myseq <- list()
  for(j in 1:length(seqfiles)){
    # define seq file
    seqfile <- read.table(seqfiles[j])
    targseq <- seqfile[seqfile$Ensembl_id == target, ]
    myseq[[j]] <- 
    
  }
  
  
  # this doesn't format correctly, might actually be better to leave as list and write as fasta
  genefile <- ldply(myseq)
  write.fasta(sequences=genefile$Sequence, names=c(genefile$Contig_id, genefile$Ensembl_id), 
              file.out=paste("genefastas/",target,".fasta", sep=""))
}


# want to save this as fasta - need to format as fasta. 
  # Use write.fasta(sequences=seqs$sequences, names=seqs$Contig_id)...something like that






