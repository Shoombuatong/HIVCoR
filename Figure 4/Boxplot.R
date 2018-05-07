
#######set directory
setwd('D:\\CRF_01AE')
#######Load package
library(randomForest)
library(kernlab)
library(e1071)
library(randomForest)
library(protr)
library(seqinr)
library(AUC)
library(ROCR)
library(RWeka)
library(caret)

######### Extract feature
x1 <- read.fasta('CCR5 AA.fasta', seqtype="AA", as.string = TRUE)
x2 <- read.fasta('CXCR4 AA.fasta', seqtype="AA", as.string = TRUE)
x1 <- x1[(sapply(x1, protcheck))]
x2 <- x2[(sapply(x2, protcheck))]
AAC1 <- t(sapply(x1, extractAAC))
AAC2 <- t(sapply(x2, extractAAC))
AAC3 <-AAC1[!duplicated(AAC1), ]
AAC4 <-AAC2[!duplicated(AAC2), ]

label =  c(rep("R5",nrow(AAC3)),rep("X4",nrow(AAC4)))
data = data.frame(rbind(AAC3,AAC4),Class = label)
#########

par( mfrow = c(3,3 ),mai=c(0.5,0.5,0.3,0.3))
colors <- c("brown1","cyan2")

boxplot(R ~ Class,data= data, main='R', col = colors,horizontal=FALSE, cex.axis=1.5,cex.main=2)
boxplot(Y ~ Class,data= data, main='Y', col = colors,horizontal=FALSE, cex.axis=1.5,cex.main=2)
boxplot(N ~ Class,data= data, main='N', col = colors,horizontal=FALSE, cex.axis=1.5,cex.main=2)
boxplot(K ~ Class,data= data, main='K', col = colors,horizontal=FALSE, cex.axis=1.5,cex.main=2)
boxplot(D ~ Class,data= data, main='D', col = colors,horizontal=FALSE, cex.axis=1.5,cex.main=2)
boxplot(S ~ Class,data= data, main='S', col = colors,horizontal=FALSE, cex.axis=1.5,cex.main=2)
boxplot(I ~ Class,data= data, main='I', col = colors,horizontal=FALSE, cex.axis=1.5,cex.main=2)
