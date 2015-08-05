library(biomaRt)
ensembl = useMart("ensembl", dataset="tguttata_gene_ensembl")
filters = listFilters(ensembl)
head(filters)
attributes = listAttributes(ensembl)
head(attributes)

test <- 'ENSTGUP00000000081'
getBM(attributes=c('ensembl_peptide_id', 'go_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=test, mart=ensembl)

biocLite("GO.db")
library("GO.db")
go_search <- getBM(attributes="go_id", filters="ensembl_peptide_id", values = test, mart = ensembl)
Term(go_search$go_id)
as.data.frame(Term(go_search$go_id))
#one of the GO terms is obsolete, use na.omit to circumvent
na.omit(as.data.frame(Term(go_search$go_id)))
for(go in go_search$go_id[1:5]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

#test with vector of peptides
test2 <- c('ENSTGUP00000012872', 'ENSTGUP00000000081', 'ENSTGUP00000010502', 'ENSTGUP00000006586', 'ENSTGUP00000003270',
           'ENSTGUP00000005098', 'ENSTGUP00000014035', 'ENSTGUP00000006602', 'ENSTGUP00000001084', 'ENSTGUP00000001948')
getBM(attributes=c('ensembl_peptide_id', 'go_id', 'hgnc_symbol'), filters='ensembl_peptide_id', values=test2, mart=ensembl)



