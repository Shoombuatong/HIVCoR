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
R5 <- read.fasta("CCR5 nlu.fasta",seqtype="DNA")
X4 <- read.fasta("CXCR4 nlu.fasta",seqtype="DNA")

codonR5  <- matrix(nrow = length(R5), ncol = 64)
codonX4  <- matrix(nrow = length(X4), ncol = 64)

for(i in 1:length(R5)){ 
codonR5[i,] = uco( R5[[i]], index = "rscu",NA.rscu = 0)
}

for(i in 1:length(X4)){ 
codonX4[i,] = uco( X4[[i]], index = "rscu",NA.rscu = 0)
}

codonR5_1 <-codonR5[!duplicated(codonR5), ]
codonX4_1 <-codonX4[!duplicated(codonX4), ]

label =  c(rep("R5",nrow(codonR5_1)),rep("X4",nrow(codonX4_1)))
internal = data.frame(rbind(codonR5_1,codonX4_1),Class = label)
AUC <- matrix(nrow = 5, ncol = 1)

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
