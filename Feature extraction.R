#######set directory
setwd('D:\\CRF_01AE')
#######Load package
library(protr)
library(seqinr

######### Amino acid composition (AAC)
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

######### Pseudo amino acid composition (PseAAC)
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

######### Relative synonymous codon usage frequency (RSCU)
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
