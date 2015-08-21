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

# load shared genes list
hits_list2 <- read.csv("Sharedhits_list.csv", header=T)
hits_list2 <- data.frame("Ensembl_id"=hits_list2$Ensembl_id)
head(hits_list2)

# load rbh files for each species
a_rbh <- read.table("RBH/A_rbh.out")
b_rbh <- read.table("RBH/B_rbh.out")
c_rbh <- read.table("RBH/C_rbh.out")
d_rbh <- read.table("RBH/D_rbh.out")
e_rbh <- read.table("RBH/E_rbh.out")
f_rbh <- read.table("RBH/F_rbh.out")
g_rbh <- read.table("RBH/G_rbh.out")
h_rbh <- read.table("RBH/H_rbh.out")
i_rbh <- read.table("RBH/I_rbh.out")
j_rbh <- read.table("RBH/J_rbh.out")
k_rbh <- read.table("RBH/K_rbh.out")
l_rbh <- read.table("RBH/L_rbh.out")

# load fasta files (assemblies)
a_fasta <- read.fasta("trinity_assembly/Trinity_A.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
b_fasta <- read.fasta("trinity_assembly/Trinity_B.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
c_fasta <- read.fasta("trinity_assembly/Trinity_C.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
d_fasta <- read.fasta("trinity_assembly/Trinity_D.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
e_fasta <- read.fasta("trinity_assembly/Trinity_E.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
f_fasta <- read.fasta("trinity_assembly/Trinity_F.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
g_fasta <- read.fasta("trinity_assembly/Trinity_G.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
h_fasta <- read.fasta("trinity_assembly/Trinity_H.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
i_fasta <- read.fasta("trinity_assembly/Trinity_I.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
j_fasta <- read.fasta("trinity_assembly/Trinity_J.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
k_fasta <- read.fasta("trinity_assembly/Trinity_K.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)
l_fasta <- read.fasta("trinity_assembly/Trinity_L.fasta", seqtype="DNA", as.string=T, forceDNAtolower=F)

# Function to match shared hits list to contig ids, and then use contig ids to pull out the sequences [note: the rbhhits and rhbseqs files need to be made first]
matchcontig_getseq <- function(infile_rbh, infile_fasta, samp){
  # match contig id from species rbh to shared gene list
  colnames(infile_rbh) <- c("contig_id", "Ensembl_id")
  hits2 <- join(hits_list, infile_rbh, by="Ensembl_id")
  # modify contig id so that it matches name in assembly fasta file
  replacebreak <- gsub(pattern="\\|", replacement=" ", x=hits2$contig_id)
  splitbreak <- strsplit(replacebreak, split=" ")
  contigname <- sapply(splitbreak, function(x){
    paste(x[[1]])
  })
  # add modified contig id to df
  hits3 <- data.frame("Ensembl_id"=hits2$Ensembl_id, "Contig_id_long"=hits2$contig_id, "Contig_id"=contigname)
  write.table(hits3, paste("rbhhits/",samp,"_hits3.txt", sep=""))
  # make fasta list to df
  fastadf <- ldply(.data=infile_fasta)
  # save file
  colnames(fastadf) <- c("Contig_id", "Sequence")
  # match contig id from rbh shared hits for given species to contig id in assembly file
  seqs <- join(hits3, fastadf, by="Contig_id") 
  # save file
  write.table(seqs, paste("rbhseqs/",samp,"_seqs.txt", sep=""))
}

matchcontig_getseq(a_rbh, a_fasta, "a")
matchcontig_getseq(b_rbh, b_fasta, "b")
matchcontig_getseq(c_rbh, c_fasta, "c")
matchcontig_getseq(d_rbh, d_fasta, "d")
matchcontig_getseq(e_rbh, e_fasta, "e")
matchcontig_getseq(f_rbh, f_fasta, "f")
matchcontig_getseq(g_rbh, g_fasta, "g")
matchcontig_getseq(h_rbh, h_fasta, "h")
matchcontig_getseq(i_rbh, i_fasta, "i")
matchcontig_getseq(j_rbh, j_fasta, "j")
matchcontig_getseq(k_rbh, k_fasta, "k")
matchcontig_getseq(l_rbh, l_fasta, "l")

# ---------------------------------------------------------
# Match up the sequences for each gene from the 12 species
# ---------------------------------------------------------

# load shared genes list
hits_list2 <- read.csv("Sharedhits_list.csv", header=T)
hits_list2 <- data.frame("Ensembl_id"=hits_list2$Ensembl_id)
head(hits_list2)


# TODO: can probably wrap this into the function
# for each of the 12 seq files, grab all rows that match given ensembl id from hits_list2
  # 
# want to save this as fasta - need to format as fasta. 
  # Use write.fasta(sequences=seqs$sequences, names=seqs$Contig_id)...something like that






