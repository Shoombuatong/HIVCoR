
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

customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

######### Extract feature
x1 <- read.fasta('CCR5 AA.fasta', seqtype="AA", as.string = TRUE)
x2 <- read.fasta('CXCR4 AA.fasta', seqtype="AA", as.string = TRUE)
x1 <- x1[(sapply(x1, protcheck))]
x2 <- x2[(sapply(x2, protcheck))]

col = 20+ 30
PAAC1  <- matrix(nrow = length(x1), ncol = col)
for (i in 1:length(x1)){
PAAC1[i,] = extractPAAC(x1[[i]][1], props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"),lambda = 30, w = 0.05, customprops = NULL)
}

PAAC2  <- matrix(nrow = length(x2), ncol = col)
for (i in 1:length(x2)){
PAAC2[i,] = extractPAAC(x2[[i]][1], props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"),lambda = 30, w = 0.05, customprops = NULL)
}

PAAC3 <-PAAC1[!duplicated(PAAC1), ]
PAAC4 <-PAAC2[!duplicated(PAAC2), ]

label =  c(rep("R5",nrow(PAAC3)),rep("X4",nrow(PAAC4)))
data2 = data.frame(rbind(PAAC3,PAAC4),Class = label)

Neg = subset(data2, Class == 'X4')
Pos = subset(data2, Class == 'R5')
nPos = nrow(Pos)
nNeg = nrow(Neg)

m= 100
ACCtr  <- matrix(nrow = m, ncol = 1)
SENStr  <- matrix(nrow = m, ncol = 1)
SPECtr  <- matrix(nrow = m, ncol = 1)
MCCtr  <- matrix(nrow = m, ncol = 1)
AUCtr  <- matrix(nrow = m, ncol = 1)
ACCts  <- matrix(nrow = m, ncol = 1)
SENSts  <- matrix(nrow = m, ncol = 1)
SPECts  <- matrix(nrow = m, ncol = 1)
MCCts  <- matrix(nrow = m, ncol = 1)
AUCts  <- matrix(nrow = m, ncol = 1)

for (i in 1:m){
#######  Dividing Training and Testing sets on positive and negative classes
sample1 <- c(sample(1:nPos,30))
sample2 <- c(sample(1:nNeg,30))
  train1  <- Pos[sample1,] ####Positive set for training
  train2  <- Neg[sample2,] ####Negative set for training
  test1 <-   Pos[-sample1,]    ####Positive set for testing
  test2 <-   Neg[-sample2,]    ####Negative set for testing 
  internal <- rbind(train1,train2) ####combining for internal set
  external <- rbind(test1,test2)    ####combining for external set

######### External validation over 5CV
tunegrid <- expand.grid(.mtry=c(1:5), .ntree=seq(100,500,100))
RFmodel <- train(Class~., data=internal , method=customRF, metric=c("Accuracy"), tuneGrid=tunegrid, trControl=trainControl(method = "cv", number=5))
Resultcv = RFmodel$ finalModel$ confusion [,1:2]
pred=prediction(RFmodel$ finalModel$ votes[,2],internal[,ncol(internal)])
perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUCtr[i,] = perf_AUC@y.values[[1]]

################### External validation
predcv = table(external$Class, predict(RFmodel, external))  ###### Prediction on external set
Resultext <- rbind(predcv[1], predcv[3],predcv[2], predcv[4]) ###### Reporting TN,FP,FN,TP
pred= prediction(predict(RFmodel ,external,type = "prob")[,2],external[,ncol(external)])
perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUCts[i,] = perf_AUC@y.values[[1]]
################### Performance report
data = Resultcv
	ACCtr[i,] = (data[1]+data[4])/(data[1]+data[2]+data[3]+data[4])*100
	SENStr[i,]  =  (data[1]/(data[1]+data[2]))*100
	SPECtr[i,] = (data[4])/(data[3]+data[4])*100
	MCC1      = (data[1]*data[4]) - (data[3]*data[2])
	MCC2      =  (data[4]+data[2])*(data[4]+data[3])
	MCC3      =  (data[1]+data[2])*(data[1]+data[3])
	MCC4	=  sqrt(MCC2)*sqrt(MCC3)
	MCCtr[i,]  = MCC1/MCC4
data = Resultext
	ACCts[i,] = (data[1]+data[4])/(data[1]+data[2]+data[3]+data[4])*100
	SENSts[i,]  =  (data[1]/(data[1]+data[2]))*100
	SPECts[i,] = (data[4])/(data[3]+data[4])*100
	MCC1      = (data[1]*data[4]) - (data[3]*data[2])
	MCC2      =  (data[4]+data[2])*(data[4]+data[3])
	MCC3      =  (data[1]+data[2])*(data[1]+data[3])
	MCC4	=  sqrt(MCC2)*sqrt(MCC3)
	MCCts[i,]  = MCC1/MCC4
}

result = data.frame (ACCtr,SPECtr,SENStr,MCCtr,AUCtr,ACCts,SPECts,SENSts,MCCts,AUCts)

Mean  <- matrix(nrow = 10, ncol = 1)
SD  <- matrix(nrow = 10, ncol = 1)
for (i in 1:10){
Mean[i,] = mean(result[,i])
SD[i,] =  sd(result[,i])
}

round(data.frame (Mean,SD),2)
