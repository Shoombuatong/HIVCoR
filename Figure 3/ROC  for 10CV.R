
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
library(Interpol)
library(Peptides)
library(reshape)

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

########################### RF-AAC
x1 <- read.fasta('CCR5 AA.fasta', seqtype="AA", as.string = TRUE)
x2 <- read.fasta('CXCR4 AA.fasta', seqtype="AA", as.string = TRUE)
x1 <- x1[(sapply(x1, protcheck))]
x2 <- x2[(sapply(x2, protcheck))]
AAC1 <- t(sapply(x1, extractAAC))
AAC2 <- t(sapply(x2, extractAAC))
AAC3 <-AAC1[!duplicated(AAC1), ]
AAC4 <-AAC2[!duplicated(AAC2), ]
label =  c(rep("R5",nrow(AAC3)),rep("X4",nrow(AAC4)))
internal = data.frame(rbind(AAC3,AAC4),Class = label)
tunegrid <- expand.grid(.mtry=c(1:5), .ntree=seq(100,500,100))
RFmodel <- train(Class~., data=internal , method=customRF, metric=c("Accuracy"), tuneGrid=tunegrid, trControl=trainControl(method = "cv", number=5))
######Loop for 10-fold CV
k <- 10;
folds <- cvsegments(nrow(internal), k);
true <- data.frame()
label <- data.frame()
for (fold in 1:k){
  currentFold <- folds[fold][[1]];
  RF = randomForest(Class ~ ., internal[-currentFold,], ntree= as.numeric(RFmodel$ bestTune[2]) ,mtry = as.numeric(RFmodel$ bestTune[1]),orm.votes=TRUE,keep.forest=TRUE, importance=TRUE) ## Building RF model
  true = rbind(true, predict(RF, internal[currentFold,],type="prob")[,2])
  label = rbind(label, internal[currentFold,]$Class)
  }
probcv= data.frame(melt(true)[,2])[,1]
write.csv(as.matrix(melt(label)$value), "matrix.csv", row.names=TRUE, na="")
label = read.csv("matrix.csv", header = TRUE)
labelcv = label[,2]

RF_aac = cbind(probcv,labelcv)

########################### RF-PseAAC
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
tunegrid <- expand.grid(.mtry=c(1:5), .ntree=seq(100,500,100))
RFmodel <- train(Class~., data=internal , method=customRF, metric=c("Accuracy"), tuneGrid=tunegrid, trControl=trainControl(method = "cv", number=5))

######Loop for 10-fold CV
k <- 10;
folds <- cvsegments(nrow(internal), k);
true <- data.frame()
label <- data.frame()

for (fold in 1:k){
  currentFold <- folds[fold][[1]];
  RF = randomForest(Class ~ ., internal[-currentFold,], ntree= as.numeric(RFmodel$ bestTune[2]) ,mtry = as.numeric(RFmodel$ bestTune[1]),orm.votes=TRUE,keep.forest=TRUE, importance=TRUE) ## Building RF model
  true = rbind(true, predict(RF, internal[currentFold,],type="prob")[,2])
  label = rbind(label, internal[currentFold,]$Class)
  }
probcv= data.frame(melt(true)[,2])[,1]
write.csv(as.matrix(melt(label)$value), "matrix.csv", row.names=TRUE, na="")
label = read.csv("matrix.csv", header = TRUE)
labelcv = label[,2]

RF_paac = cbind(probcv,labelcv)

########################### RF-RSCU
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
tunegrid <- expand.grid(.mtry=c(1:5), .ntree=seq(100,500,100))
RFmodel <- train(Class~., data=internal , method=customRF, metric=c("Accuracy"), tuneGrid=tunegrid, trControl=trainControl(method = "cv", number=5))

######Loop for 10-fold CV
k <- 10;
folds <- cvsegments(nrow(internal), k);
true <- data.frame()
label <- data.frame()

for (fold in 1:k){
  currentFold <- folds[fold][[1]];
  RF = randomForest(Class ~ ., internal[-currentFold,], ntree= as.numeric(RFmodel$ bestTune[2]) ,mtry = as.numeric(RFmodel$ bestTune[1]),orm.votes=TRUE,keep.forest=TRUE, importance=TRUE) ## Building RF model
  true = rbind(true, predict(RF, internal[currentFold,],type="prob")[,2])
  label = rbind(label, internal[currentFold,]$Class)
  }
probcv= data.frame(melt(true)[,2])[,1]
write.csv(as.matrix(melt(label)$value), "matrix.csv", row.names=TRUE, na="")
label = read.csv("matrix.csv", header = TRUE)
labelcv = label[,2]

RF_rscu = cbind(probcv,labelcv)



