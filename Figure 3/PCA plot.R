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

#####################################################
res.pca <- PCA(data[,-n])

write.csv(res.pca$eig, "Eig.csv", row.names=FALSE, na="")
write.csv(res.pca$var$coord, "Loading.csv", row.names=FALSE, na="")
write.csv(res.pca$ind$coord, "Score.csv", row.names=FALSE, na="")

score_active1 = read.csv("Score.csv", header = TRUE)
loading_active1 = read.csv("Loading.csv", header = TRUE)

#####################################################
par( mfrow = c(1,2 ),mai=c(0.9,0.9,0.9,0.9))

Myred <- rgb(t(col2rgb("red")), alpha=150, maxColorValue=300)
Myblue <- rgb(t(col2rgb("blue")), alpha=150, maxColorValue=300)

plot(Dim.1 ~ Dim.2, score_active1, main= NULL,xlab = 'PC1(16.16%)' ,ylab = 'PC2(10.52%)', type = "n",pch=16, lwd =5, xlim=c(-6, 6), ylim=c(-6, 6),
cex.axis=1.5,cex.lab=1.2)

data1 = subset(score_active1, Class == 'R5')
points(data1[,1], data1[,2], col = Myred,pch=16, lwd =2,cex = 1.5)  
data1 = subset(score_active1, Class == 'X4')
points(data1[,1], data1[,2], col = Myblue,pch=16, lwd =2,cex = 1.5)  


plot(loading_active1[,1], loading_active1[,2], main= NULL, xlab= "PC1(16.16%)",ylab= "PC2(10.52%)", col= "forestgreen",pch=19, xlim=c(-1, 1), ylim=c(-1, 1),cex.axis=1.5,cex.lab=1.2, lwd =3)
text(loading_active1[,1], loading_active1[,2], labels= loading_active1[,3], cex= 1.2,pos=3)
abline(h=0)
abline(v=0)
