setwd("Data/")

# LOAD ======
set.seed(5099)
# load(file = "Data/.RData")

allecc <- read.csv("All data All cohorts.csv", header = T)
                   # , check.names=F)
allprox <- read.csv("Proximity All Cohorts.csv", header = T)
                    # , check.names=F)
names(allprox)[1] <- paste("CaseNo")
data <- merge(allecc, allprox, by = "CaseNo", sort = F)

# remove duplicate columns
data <- data[,-c(156:157)]

# set DiseaseSpecificDeath, Sex, Age, pT, Site, Diff as factor
data[,c(64:70)] <- lapply(data[,c(64:70)], as.factor)
data$set <- ""
data$set[1:113] <- "Train"
data$set[114:169] <- "Validation(EDI)"
data$set[170:230] <- "Validation(JAP)"

# training and validation dataset
data1 <- data
# drop normalised columns
data1 <- data[, -grep("norm", colnames(data))]

# replace NA with mode, i.e. 1 for Diff
data1$Diff[228] <-1

train <- data1[1:113, -c(1,126)]
v1 <- data1[114:169, -c(1,126)]
v2 <- data1[170:230, -c(1,126)]




# STANDARDIZATION ======
# standardize except DiseaseSpecificSurvival and categorical var
data1[,-c(1, 63:70,126)] <- scale(data1[,-c(1, 63:70,126)], center = T, scale = T)

# training and validation dataset
trains <- data1[1:113, -c(1,126)]
v1s <- data1[114:169, -c(1,126)]
v2s <- data1[170:230, -c(1,126)]

# extract design matrix 
chop <- trains[,-c(62:63)]
chop.v1 <- v1s[,-c(62:63)]
chop.v2 <- v2s[,-c(62:63)]

# dummy for glmnet
xfactorv1 <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=chop.v1)[, -1]
chop.v1 <- as.matrix(data.frame(chop.v1[-c(62:67)], xfactorv1))
chop.v1 <- data.frame(chop.v1)

xfactorv2 <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=chop.v2)[,-1]
chop.v2 <- as.matrix(data.frame(chop.v2[-c(62:67)], xfactorv2))
chop.v2 <- data.frame(chop.v2)
