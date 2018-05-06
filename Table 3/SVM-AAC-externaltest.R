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

label =  c(rep("R5",150),rep("X4",40))
data2 = data.frame(rbind(AAC3,AAC4),Class = label)

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

######### Optimized parameter
SVMpara = tune(svm, Class ~ ., data =internal, ranges =list(gamma = 2^(-8:8), cost = 2^(-8:8)),
tunecontrol = tune.control(sampling = "fix"))

################### 10-fold CV
k <- 10;
Resultcv <- 0;
folds <- cvsegments(nrow(internal), k);
for (fold in 1:k){
    currentFold <- folds[fold][[1]];
    RF <- ksvm(Class~.,data=internal[-currentFold,],kernel="rbfdot", cost = as.numeric(SVMpara$ best.parameters[1]), gamma = as.numeric(SVMpara$ best.parameters[1]),prob.model=TRUE)
    pred = predict(RF, internal[currentFold,])
    Resultcv <- Resultcv + table(true=internal[currentFold,]$Class, pred=pred);
pred2=prediction(predict(RF, internal[currentFold,], type="probabilities")[,2],internal[currentFold,]$Class)
perf_AUC=performance(pred2,"auc") #Calculate the AUC value
AUC[fold,] = perf_AUC@y.values[[1]]
}
AUCtr[i,] = mean(AUC)
################### External validation
RF = ksvm(Class~.,data=internal,kernel="rbfdot", cost = as.numeric(SVMpara$ best.parameters[1]), gamma = as.numeric(SVMpara$ best.parameters[1]),prob.model=TRUE)
predcv = table(external$Class, predict(RF, external))  ###### Prediction on external set
Resultext <- rbind(predcv[1], predcv[3],predcv[2], predcv[4]) ###### Reporting TN,FP,FN,TP
pred2=prediction(predict(RF, external, type="probabilities")[,2],external$Class)
perf_AUC=performance(pred2,"auc") #Calculate the AUC value
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
