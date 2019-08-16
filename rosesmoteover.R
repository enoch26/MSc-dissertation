# setwd("Data/")

# Output 
# overchop  - oversample standardized design matrix with dummy variables for glmnet(smote, rose or ovunsample)
# data4     - oversample standardized data
# fold      - 10 fold cv foldid
# LOAD ======
# load(file = "Data/.RData")
# set seed
set.seed(5099)


# load("oversample.RData") 

library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
# stopCluster(cl)

# library ====
library(ROSE)
library(DMwR)
library(caret)
library(MASS)
library(dplyr)
library(ggplot2)
library(kernlab)
library(klaR)
library(biotools)
library(survival)
library(survminer)
library(gbm)
library(randomForest)
library(randomForestSRC)
library(ggRandomForests)
library(naivebayes)
library(Hmisc)
library(glmnet)

# data input =====
source("data input.R")

# oversampling ======
# train data with ovun.sample, ROSE and SMOTE
trains$DiseaseSpecificDeath <- as.factor(trains$DiseaseSpecificDeath)

# data.ovun <- ovun.sample(DiseaseSpecificDeath ~ ., 
#                                   data = trains, method = "both", 
#                                   p=0.5, N = 1000, seed=5099)$data
# table(data.ovun$DiseaseSpecificDeath)
# 
# data.rose <- ROSE(DiseaseSpecificDeath ~ ., data = trains, N=1000, seed = 5099)$data
# table(data.rose$DiseaseSpecificDeath)

data.smote <- SMOTE(DiseaseSpecificDeath ~., data = trains, perc.over = 3200, perc.under = 105, k = 5)
as.data.frame(table(data.smote$DiseaseSpecificDeath))

# choose ovun.sample, rose or smote =========
# data4 <- data.ovun/data.rose/data.smote
data4 <- data.smote
overchop <- data.smote[,-c(62:63)]

# create dummy for glmnet
xfactor <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=overchop)[, -1]
# combine
overchop <- as.matrix(data.frame(overchop[-c(61:66)], xfactor))
overchop <- data.frame(overchop)
names(overchop)

# for reproducible 10 fold cv.glmnet
ind <- sample(1:nrow(overchop),nrow(overchop))
data.sca <- data.sca[ind,]
fold <- rep(1:10,length.out=nrow(overchop))

save.image(file='rosesmoteover.RData')
# sessioninfo =========
# sessionInfo()
