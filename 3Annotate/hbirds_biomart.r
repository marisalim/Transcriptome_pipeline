# Marisa Lim (c)2015
MLwd = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files"
setwd(MLwd)

# Load libraries
require(plyr)
library(seqinr)

# -------------------------------------
# find RBHs shared between all species
# -------------------------------------
version2 <- function(){
  # Create infile lists <- these are the outputs from theblastparser.py
  rbhfolder = "RBH/"
  rbhfiles <- list.files(paste(rbhfolder,sep="/"),pattern=".out",full.name=T,recursive=T)
  
  complist <- list()
  for(i in 1:length(rbhfiles)){
    # load infile_rbh
    infile_rbh <- read.table(rbhfiles[i])
    colnames(infile_rbh) <- c("contig_id", "Ensembl_id", "Query_start", "Query_end", "Target_start", "Target_end", "evalue", "Bit_score")
    complist[[i]] <- as.vector(infile_rbh[,2])
  }
  hits <- Reduce(intersect, complist)
  #length(hits)
  #head(hits)
  hitsvec <- as.vector(hits)
  #length(hitsvec)
}
version2()

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
# =====================
grep("go",attributes$description, ignore.case=T, value=T)

test <- 'ENSTGUP00000000081'
getBM(attributes=c('ensembl_peptide_id', 'go_id','goslim_goa_description','name_1006','definition_1006', 
                   'go_linkage_type', 'namespace_1003', 'goslim_goa_accession',
                   'hgnc_symbol'), filters='ensembl_peptide_id', values=test, mart=ensembl)

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

# ===========================
# now biomaRt with real data:
# ===========================
hits_go <- getBM(attributes=c('ensembl_peptide_id', 'go_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=hitsvec, mart=ensembl)
head(hits_go)
dim(hits_go)
class(hits_go)

goterms <- as.data.frame(Term(hits_go$go_id))
head(goterms)
class(goterms)

gotermstogether <- data.frame("ensembl_peptide_id"=hits_go$ensembl_peptide_id, "hgnc_symbol"=hits_go$hgnc_symbol, "go_id"=hits_go$go_id,
                              "go_description"=goterms$`Term(hits_go$go_id)`)
# gotermstogether[gotermstogether$hgnc_symbol == "TPP1", ]
dat_tab <- table(gotermstogether[, 4]) 
write.table(dat_tab, "Fulldatasetgenes_gotermstable.txt")
barplot(dat_tab)

hits_go2 <- getBM(attributes=c('ensembl_peptide_id', 'hgnc_symbol', 'description'), filters='ensembl_peptide_id', values=hitsvec, mart=ensembl)
head(hits_go2)
dim(hits_go2)

write.csv(hits_go2, "Sharedhits_gene_descript.csv")

hits_list <- data.frame(Ensembl_id = hits_go2[,1])
write.csv(hits_list, "Sharedhits_list.csv")

# -------------------------------------------------------------------------
# Match ensembl id of shared hits to contig ID (need this to find sequence
# -------------------------------------------------------------------------
# Function to match shared hits list to contig ids
matchcontig <- function(){
  # make output folder
  outfile = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files/rbhhits"
  dir.create(outfile)
  
  # load shared genes list
  hits_list2 <- read.csv("Sharedhits_list.csv", header=T)
  hits_list2 <- data.frame("Ensembl_id"=hits_list2$Ensembl_id)
  #head(hits_list2)
  
  # Create infile lists
  rbhfolder = "RBH/"
  rbhfiles <- list.files(paste(rbhfolder,sep="/"),pattern=".out",full.name=T,recursive=T)
  
  # Create sample id vector
  sampleids <- c("a_Phart", "b_Cmuls", "c_Pmala", "d_Acast", "e_Ccoru", "f_Ccoel", "g_Cviol", "h_Aamaz", "i_Mphoe", "j_Aviri", "k_Amela", "l_Pgiga")
  
  # Match shared hits to species rbh
  for(i in 1:length(rbhfiles)){
    # load infile_rbh
    infile_rbh <- read.table(rbhfiles[i])
    colnames(infile_rbh) <- c("contig_id", "Ensembl_id", "Query_start", "Query_end", "Target_start", "Target_end", "evalue", "Bit_score")
    
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
    hits3 <- data.frame("Ensembl_id"=hits2$Ensembl_id, "Contig_id_long"=hits2$contig_id, "Contig_id"=contigname, "Sample_id"=hits2$Sample_id, 
                        "Query_start"=hits2$Query_start, "Query_end"=hits2$Query_end, "Target_start"=hits2$Target_start, 
                        "Target_end"=hits2$Target_end, "evalue"=hits2$evalue, "Bit_score"=hits2$Bit_score)
    write.table(hits3, paste("rbhhits/",samp,"_hits3.txt", sep=""))
  }
}
matchcontig()

# --------------------------------------
# Get sequences from matched contig ids
# --------------------------------------
# Function that uses contig ids to pull out the sequences
getmyseqs <- function(){
  # make output folder
  outfile = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files/rbhseqs"
  dir.create(outfile)
  
  # load hits files for each species
  hits3folder <- "rbhhits/"
  hits3files <- list.files(paste(hits3folder, sep="/"), pattern=".txt", full.name=T, recursive=T)
  
  # Create infile lists
  fastafolder = "trinity_assembly/"
  fastafiles <- list.files(paste(fastafolder,sep="/"),pattern=".fasta",full.name=T,recursive=T)
  
  # Create sample id vector
  sampleids <- c("a_Phart", "b_Cmuls", "c_Pmala", "d_Acast", "e_Ccoru", "f_Ccoel", "g_Cviol", "h_Aamaz", "i_Mphoe", "j_Aviri", "k_Amela", "l_Pgiga")
  
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
getmyseqs()

# ---------------------------------------------------------
# Match up the sequences for each gene from the 12 species
# ---------------------------------------------------------

# for each of the 12 seq files, grab all rows that match given ensembl id from hits_list2
makegenefiles <- function(){
  # make output folder
  outfile = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files/genefastas"
  dir.create(outfile)
  
  # load shared genes list
  hits_list2 <- read.csv("Sharedhits_list.csv", header=T)
  hits_list2 <- data.frame("Ensembl_id"=hits_list2$Ensembl_id)
  #head(hits_list2)
  
  # load seqs files for each sample
  seqfolder = "rbhseqs/"
  seqfiles <- list.files(paste(seqfolder,sep="/"),pattern=".txt",full.name=T,recursive=T)
  
  for(i in 1:nrow(hits_list2)){
    # define ensembl id to find
    target <- hits_list2[i,]
    myseq <- list()
    
    # find target in seq files
    for(j in 1:length(seqfiles)){
      # define seq file
      seqfile <- read.table(seqfiles[j])
      myseq[[j]] <- seqfile[seqfile$Ensembl_id == target, ]
    }
    genefile <- ldply(myseq)
    genelist <- as.list(genefile$Sequence)
    
    # save the sequences for each gene
    write.fasta(sequences=genelist, 
                names=paste(genefile$Sample_id,".",genefile$Ensembl_id,".",genefile$Contig_id,sep=""), 
                nbchar=80, file.out=paste("genefastas/",target,".fasta", sep=""))
    # save the blast info for each gene, rearrange columns to be in same order as original blast output
    genefile2 <- genefile[,c(3,1,4:10)]
    write.table(genefile2, file=paste("genefastas/",target,"_blastout.txt", sep=""))
  }
}
makegenefiles()







