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
library(pls)
internal = read.csv("number training set R5 X4 SVMLMT.csv", header = TRUE) 
AUC <- matrix(nrow = 10, ncol = 1)

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
AUCtr= mean(AUC)

data = Resultcv
	ACCtr = (data[1]+data[4])/(data[1]+data[2]+data[3]+data[4])*100
	SENStr  =  (data[1]/(data[1]+data[2]))*100
	SPECtr = (data[4])/(data[3]+data[4])*100
	MCC1      = (data[1]*data[4]) - (data[3]*data[2])
	MCC2      =  (data[4]+data[2])*(data[4]+data[3])
	MCC3      =  (data[1]+data[2])*(data[1]+data[3])
	MCC4	=  sqrt(MCC2)*sqrt(MCC3)
	MCCtr  = MCC1/MCC4

round(data.frame (ACCtr,SPECtr,SENStr,MCCtr,AUCtr),2)
