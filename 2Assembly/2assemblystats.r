# Code to look at the contig length distributions for each Trinity assembly
# Note: I had to redo the trinity assemblies, but I've left in the old contig distributions for comparison for now
# Marisa Lim (c)2015

# set working directory
wd = "C:/Users/mcwlim/Dropbox/Marisacompfiles/Transcriptome files/"
setwd(wd)

# Load contig length input data
contigA <- read.table("extractidsA.txt")
contigA_redo <- read.table("extract_contiglen_A.txt")
contigB <- read.table("extractidsB.txt")
contigB_redo <- read.table("extract_contiglen_B.txt")
contigC <- read.table("extractidsC.txt")
contigC_redo <- read.table("extract_contiglen_C.txt")
contigD <- read.table("extractidsD.txt")
contigD_redo <- read.table("extract_contiglen_D.txt")
contigE <- read.table("extractidsE.txt")
contigE_redo <- read.table("extract_contiglen_E.txt")
contigF <- read.table("extractidsF.txt")
contigF_redo <- read.table("extract_contiglen_F.txt")
contigG <- read.table("extractidsG.txt")
contigH <- read.table("extractidsH.txt")
contigH_redo <- read.table("extract_contiglen_H.txt")
contigI <- read.table("extractidsI.txt")
contigI_redo <- read.table("extract_contiglen_I.txt")
contigJ <- read.table("extractidsJ.txt")
contigJ_redo <- read.table("extract_contiglen_J.txt")
contigK <- read.table("extractidsK.txt")
contigL <- read.table("extractidsL.txt")
contigL_redo <- read.table("extract_contiglen_L.txt")

contiglen <- sort(contigA$V1) # if you want to sort it..
hist(as.numeric(contigA$V1))

# plot distribution for each sample
jpeg("CGML001A-L_contigdistr.jpg", height=10, width=15, units="in", res=600)
par(mfrow=c(4,6))
hist(as.numeric(contigA$V1), breaks=50, col="turquoise", xlab="Contig length", 
     main="CGML001A")
hist(as.numeric(contigA_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001A_redo')
hist(as.numeric(contigB$V1), breaks=50, col="turquoise", xlab="Contig length", 
     main="CGML001B")
hist(as.numeric(contigB_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001B_redo')
hist(as.numeric(contigC$V1), breaks=50, col="turquoise", xlab="Contig length", 
     main="CGML001C")
hist(as.numeric(contigC_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001C_redo')
hist(as.numeric(contigD$V1), breaks=50, col="turquoise", xlab="Contig length", 
     main="CGML001D")
hist(as.numeric(contigD_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001D_redo')
hist(as.numeric(contigE$V1), breaks=50, col="turquoise", xlab="Contig length", 
     main="CGML001E")
hist(as.numeric(contigE_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001E_redo')
hist(as.numeric(contigF$V1), breaks=50, col="turquoise", xlab="Contig length", 
     main="CGML001F")
hist(as.numeric(contigF_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001F_redo')
hist(as.numeric(contigH$V1), breaks=50, col="turquoise", xlab="Contig length",
     main="CGML001H")
hist(as.numeric(contigH_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001H_redo')
hist(as.numeric(contigI$V1), breaks=50, col="turquoise", xlab="Contig length",
     main="CGML001I")
hist(as.numeric(contigI_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001I_redo')
hist(as.numeric(contigJ$V1), breaks=50, col="turquoise", xlab="Contig length",
     main="CGML001J")
hist(as.numeric(contigJ_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001J_redo')
hist(as.numeric(contigL$V1), breaks=50, col="turquoise", xlab="Contig length",
     main="CGML001L")
hist(as.numeric(contigL_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001L_redo')
hist(as.numeric(contigG$V1), breaks=50, col="turquoise", xlab="Contig length",
     main="CGML001G")
hist(as.numeric(contigK$V1), breaks=50, col="turquoise", xlab="Contig length", 
     main="CGML001K")
dev.off()

jpeg("CGML001A-L_contiglen_4April16.jpg", height=10, width=10, units="in", res=600)
par(mfrow=c(3,4))
hist(as.numeric(contigA_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001A_redo')
hist(as.numeric(contigB_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001B_redo')
hist(as.numeric(contigC_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001C_redo')
hist(as.numeric(contigD_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001D_redo')
hist(as.numeric(contigE_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001E_redo')
hist(as.numeric(contigF_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001F_redo')
hist(as.numeric(contigH_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001H_redo')
hist(as.numeric(contigI_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001I_redo')
hist(as.numeric(contigJ_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001J_redo')
hist(as.numeric(contigL_redo$V1), breaks=50, col="red", xlab="Contig length", main='CGML001L_redo')
dev.off()
