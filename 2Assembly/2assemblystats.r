# Code to look at the contig length distributions for each Trinity assembly
# Marisa Lim (c)2015

# set working directory
wd = "C:/Users/mcwlim/Desktop/StonyBrook/GrahamLab/Dissertation idea materials/THESIS PROJECTS/Hummingbird_Genomics/Transcriptome pipeline files/"
setwd(wd)

# Load contig length input data
contigA <- read.table("extractidsA.txt")
contigB <- read.table("extractidsB.txt")
contigC <- read.table("extractidsC.txt")
contigD <- read.table("extractidsD.txt")
contigE <- read.table("extractidsE.txt")
contigF <- read.table("extractidsF.txt")
contigG <- read.table("extractidsG.txt")
contigH <- read.table("extractidsH.txt")
contigI <- read.table("extractidsI.txt")
contigJ <- read.table("extractidsJ.txt")
contigK <- read.table("extractidsK.txt")
contigL <- read.table("extractidsL.txt")

contiglen <- sort(contigA$V1) # if you want to sort it..
hist(as.numeric(contigA$V1))

# plot distribution for each sample
jpeg("CGML001A-L_contigdistr.jpg", height=10, width=15, units="in", res=600)
par(mfrow=c(3,4))
hist(as.numeric(contigA$V1), breaks=50, col="tomato", xlab="Contig length", 
     main="CGML001A")
hist(as.numeric(contigB$V1), breaks=50, col="indianred", xlab="Contig length", 
     main="CGML001B*")
hist(as.numeric(contigC$V1), breaks=50, col="violetred", xlab="Contig length", 
     main="CGML001C")
hist(as.numeric(contigD$V1), breaks=50, col="salmon2", xlab="Contig length", 
     main="CGML001D")
hist(as.numeric(contigE$V1), breaks=50, col="orange", xlab="Contig length", 
     main="CGML001E*")
hist(as.numeric(contigF$V1), breaks=50, col="gold", xlab="Contig length", 
     main="CGML001F")
hist(as.numeric(contigG$V1), breaks=50, col="yellowgreen", xlab="Contig length",
     main="CGML001G")
hist(as.numeric(contigH$V1), breaks=50, col="darkgreen", xlab="Contig length",
     main="CGML001H")
hist(as.numeric(contigI$V1), breaks=50, col="turquoise", xlab="Contig length",
     main="CGML001I")
hist(as.numeric(contigJ$V1), breaks=50, col="navy", xlab="Contig length",
     main="CGML001J*")
hist(as.numeric(contigK$V1), breaks=50, col="plum4", xlab="Contig length", 
     main="CGML001K")
hist(as.numeric(contigL$V1), breaks=50, col="purple4", xlab="Contig length",
     main="CGML001L")
dev.off()

