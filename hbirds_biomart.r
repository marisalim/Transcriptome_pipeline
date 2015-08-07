#Marisa Lim (c)2015
MLwd = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files"
setwd(MLwd)

# find RBHs shared between all species
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

#Use biomaRt to find gene names and GO info for these hits
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
