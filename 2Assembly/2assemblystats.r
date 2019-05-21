# Code to look at the contig length distributions for each Trinity assembly
# Marisa Lim (c)2015

# set working directory
wd = "[set your working directory]"
setwd(wd)

# Load contig length input data
contigA_redo <- read.table("extract_contiglen_A.txt")
contigB_redo <- read.table("extract_contiglen_B.txt")
contigC_redo <- read.table("extract_contiglen_C.txt")
contigD_redo <- read.table("extract_contiglen_D.txt")
contigE_redo <- read.table("extract_contiglen_E.txt")
contigF_redo <- read.table("extract_contiglen_F.txt")
contigG_redo <- read.table("extract_contiglen_G.txt")
contigH_redo <- read.table("extract_contiglen_H.txt")
contigI_redo <- read.table("extract_contiglen_I.txt")
contigJ_redo <- read.table("extract_contiglen_J.txt")
contigK_redo <- read.table("extract_contiglen_K.txt")
contigL_redo <- read.table("extract_contiglen_L.txt")

contiglen <- sort(contigA$V1) # if you want to sort it..
hist(as.numeric(contigA$V1))

# plot distribution for each sample
jpeg("CGML001A-L_contigdistr.jpg", height=10, width=15, units="in", res=600)
par(mfrow=c(4,3))
hist(as.numeric(contigA_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001A_redo')
hist(as.numeric(contigB_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001B_redo')
hist(as.numeric(contigC_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001C_redo')
hist(as.numeric(contigD_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001D_redo')
hist(as.numeric(contigE_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001E_redo')
hist(as.numeric(contigF_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001F_redo')
hist(as.numeric(contigG_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001G_redo')
hist(as.numeric(contigH_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001H_redo')
hist(as.numeric(contigI_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001I_redo')
hist(as.numeric(contigJ_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001J_redo')
hist(as.numeric(contigK_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001K_redo')
hist(as.numeric(contigL_redo$V1), breaks=50, col="turquoise", xlab="Contig length", main='CGML001L_redo')
dev.off()
