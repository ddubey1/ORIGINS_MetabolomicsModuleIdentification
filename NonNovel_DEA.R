#----------------- Installations -----------------------------------------------

# installing limma

#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("limma")

# installing edgeR

# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("edgeR")

# ------------------------------------------------------------------------------



# -------------------- Differential Expression Analysis ------------------------

# load libraries
library(limma)
library(edgeR)
limmaUsersGuide()

# set working directory [CHANGE DEPENDING ON USER]
getwd()
setwd("/Users/dhruv/Desktop/Summer 2022/AH Research/ORIGINS Multi-Omics : Metabolomics")

# read in raw expression data
exprs_nn <- read.csv("Exprs_NonNovel.csv")
exprs_nn <- exprs_nn[1:909,] # first 909 rows of csv are non-novel

# info
dim(exprs_nn) # 909 by 61
names(exprs_nn)

# 70% 0s filtering
zeroes = exprs_nn[,-1] == 0
zeroes_sum = rowSums(zeroes)
idx_keep = which(zeroes_sum < 0.7 * ncol(zeroes))
length(idx_keep)
exprs_nn = exprs_nn[idx_keep,] # now 907 by 61 (2 metabolites removed)

# log2 transform expression data
exprs_nn[,-1] = log(exprs_nn[,-1] + 1, base = 2)
dim(exprs_nn)

# convert expression data to matrix form
exprs_matrix = as.matrix(exprs_nn, rownames = TRUE)
rownames(exprs_matrix) = exprs_matrix[,1]
exprs_matrix = exprs_matrix[,-1]

# make numeric
exprs_matrix_num <- matrix(as.numeric(exprs_matrix), ncol = ncol(exprs_matrix))
rownames(exprs_matrix_num) = rownames(exprs_matrix)
colnames(exprs_matrix_num) = colnames(exprs_matrix)

# quantile normalization
norm_exprs_matrix = normalizeBetweenArrays(exprs_matrix_num)

# create group
sample_names = colnames(norm_exprs_matrix)
sample_names
dtype <- substr(sample_names, 1, 3) 
time <- substr(sample_names, nchar(sample_names) - 4, nchar(sample_names) - 1)
for(i in 1:length(time)) {
  if(time[i] == "GTT1"){
    time[i] = "OGTT"
  }
  if(time[i] == "ase1"){
    time[i] = "Base"
  }
}

time
dtype

group <- interaction(dtype, time)
group

# -----------------------------------------------------------------------------
# read in clinical data (for design matrix) 
clinical = read.csv("clinical_tofix.csv")

# make male 1 and females 0
rownames(clinical) = clinical[,1]
clinical = clinical[,-1]
gender_row = clinical[1,]
for(i in 1:length(gender_row)) {
  if(gender_row[i]=="F") {
    gender_row[i] = 0
  }
  else {
    gender_row[i] = 1
  }
}
clinical[1,] = gender_row

# clinical is created
clinical

gender = as.numeric(clinical[1,])
age = as.numeric(clinical[2,])
BMI = as.numeric(clinical[3,])


# create design matrix
design <- model.matrix(~ 0 + group + gender + age + BMI)
# design <- model.matrix(~ 0 + group)
design

# create contrast matrix
ctrx <- makeContrasts(Pre.vs.NGT_Base = groupPre.Base - groupNGT.Base, T2D.vs.NGT_Base = groupT2D.Base - groupNGT.Base,
                      PDM.vs.NGT_Change = (groupPre.OGTT-groupPre.Base)-(groupNGT.OGTT-groupNGT.Base),
                      T2D.vs.NGT_Change = (groupT2D.OGTT-groupT2D.Base)-(groupNGT.OGTT-groupNGT.Base), levels = design)
ctrx

# perform linear modeling
fit1 <- lmFit(norm_exprs_matrix, design = design)
fit2 <- contrasts.fit(fit1, contrast = ctrx)
fit3 <- eBayes(fit2)

PDM.NGT_Base <- topTable(fit3, coef = 1, p.value = 0.2, lfc = log2(1.5))
T2D.NGT_Base <- topTable(fit3, coef = 2, p.value = 0.2, lfc = log2(1.5))
PDM.NGT_Change <- topTable(fit3, coef = 3, p.value = 0.2, lfc = log2(1.5))
T2D.NGT_Change <- topTable(fit3, coef = 4, p.value = 0.2, lfc = log2(1.5))

?topTable




