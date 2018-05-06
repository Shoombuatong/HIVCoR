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
AUC <- matrix(nrow = 5, ncol = 1)

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
internal = data.frame(rbind(PAAC3,PAAC4),Class = label)
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
