#Marisa Lim (c)2015
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
#Use biomaRt to find gene names and GO info for these hits
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
# -------------------------------------------------------------------------
require(plyr)
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

matchcontig <- function(infile, samp){
  colnames(infile) <- c("contig_id", "Ensembl_id")
  hits2 <- join(hits_list, infile, by="Ensembl_id")
  replacebreak <- gsub(pattern="\\|", replacement=" ", x=hits2$contig_id)
  splitbreak <- strsplit(replacebreak, split=" ")
  contigname <- sapply(splitbreak, function(x){
    paste(x[[1]])
  })
  hits3 <- data.frame("Ensembl_id"=hits2$Ensembl_id, "Contig_id_long"=hits2$contig_id, "Contig_id"=contigname)
  write.table(hits3, paste(samp,"_hits3.txt", sep=""))
}
matchcontig(a_rbh, "a")
matchcontig(b_rbh, "b")
matchcontig(c_rbh, "c")
matchcontig(d_rbh, "d")
matchcontig(e_rbh, "e")
matchcontig(f_rbh, "f")
matchcontig(g_rbh, "g")
matchcontig(h_rbh, "h")
matchcontig(i_rbh, "i")
matchcontig(j_rbh, "j")
matchcontig(k_rbh, "k")
matchcontig(l_rbh, "l")

# --------------------------------
# Match contig id and get sequence
# --------------------------------
library(seqinr)
a_fasta <- read.fasta("trinity_assembly/Trinity_A.fasta", seqtype="DNA", as.string = T)
a_fasta[1]
a_seqs <- list()
for(i in 1:length(a_hits2)){
  if(a_hits2$contig_id[i] == a_fasta[i]){
    a_seqs[[i]] <- a_fasta[i]
  } else{
    i <- i + 1
  }
}



