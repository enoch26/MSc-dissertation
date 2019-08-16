#The function linreg performs a bootstrap on glmnet with lambda.min
#Input Arguments:
#       B           - # of iterations
#       X           - design matrix containing predictor variables
#       time        - vector containing survial time
#       event       - vector containing binary event
#       a           - alpha for glmnet
#       fold        - # of cv fold
#       oversample  - T: oversample n=1000 and p=0.5 via ovun.sample, package ROSE, 
#                     F: no oversample
# Outcome:
#       probzero    - probability of zero for each predictor var
#       beta        - matrix of beta in each iteration at lambda.min

source("data input.R")
save.image(file="glmnetboot.RData")
load("glmnetboot.RData")

library(ROSE)
library(glmnet)
library(survival)

# change to dummy for glmnet
xfactor <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=chop)[, -1]
chop <- as.matrix(data.frame(chop[-c(62:67)], xfactor))

trains$DiseaseSpecificDeath <- as.numeric(trains$DiseaseSpecificDeath)

glmnetboot <- function(B, X, time, event, a=1, fold, oversample=T) {
  
  beta <- matrix(0, nrow = ncol(X), ncol = B)
  nnz <- list()
  if (oversample==T) {
    data <- data.frame(X, time, event)
    for(i in 1:B){
      # oversample
      X <- ovun.sample(event ~ ., 
                       data = data, method = "both", 
                       p=0.5, N = 1000)$data
      y <- Surv(X$time, X$event)
      X <- subset(X, select = -c(time,event) )
      X <- as.matrix(X)
      # glmnet model
      glmnet_model <- glmnet(X, y, family = "cox", alpha = a)
      cv.fit <- cv.glmnet(X, y, alpha = a,
                          family="cox", nfold = fold, keep = F,
                          parallel = T)
      lambda <- cv.fit$lambda.min
      beta[,i] <- as.matrix(coef(glmnet_model, s = lambda))
    }
    rownames(beta) <- colnames(X)
    beta <- t(beta)
    nnz <- apply(beta == 0, 2, sum)/B
    return(list("probzero" = nnz, "beta" = beta))
  }
  else{
    y <- Surv(time, event)
    glmnet_model <- glmnet(X, y, family = "cox", alpha = a)
    for(i in 1:B){
      cv.fit <- cv.glmnet(X, y, alpha = a,
                          family="cox", nfold = fold, keep = F,
                          parallel = T)
      lambda <- cv.fit$lambda.min
      beta[,i] <- as.matrix(coef(glmnet_model, s = lambda))
    }
    rownames(beta) <- colnames(X)
    beta <- t(beta)
    nnz <- apply(beta == 0, 2, sum)/B
    return(list("probzero" = nnz, "beta" = beta))
  }
}
#  n=100
r <- glmnetboot(B=100, X=chop, 
                 time = trains$DiseaseSpecificSurvival, 
                 event = trains$DiseaseSpecificDeath, 
                 a=0, fold=10, oversample = T)

l <- glmnetboot(B=100, X=chop, 
                 time = trains$DiseaseSpecificSurvival, 
                 event = trains$DiseaseSpecificDeath, 
                 a=1, fold=10, oversample = T)

e <- glmnetboot(B=100, X=chop, 
                 time = trains$DiseaseSpecificSurvival, 
                 event = trains$DiseaseSpecificDeath, 
                 a=0.5, fold=10, oversample = T)

write.csv(r$beta, file="ridgebeta.csv")
write.csv(l$beta, file="lassobeta.csv")
write.csv(e$beta, file="enetbeta.csv")

par(mar = c(7, 4.1, 2, 2.1))
boxplot(r$beta, las=2,cex.axis=0.3)

boxplot(l$beta, las=2,cex.axis=0.3)
barplot(l$probzero, las=2, cex.names=0.3)

boxplot(l$beta, las=2,cex.axis=0.3)
barplot(l$probzero, las=2, cex.names=0.3)

boxplot(e$beta, las=2,cex.axis=0.3)
barplot(e$probzero, las=2, cex.names=0.3)

# n=1000
r2 <- glmnetboot(B=1000, X=chop, 
                time = trains$DiseaseSpecificSurvival, 
                event = trains$DiseaseSpecificDeath, 
                a=0, fold=10, oversample = T)

l2 <- glmnetboot(B=1000, X=chop, 
                time = trains$DiseaseSpecificSurvival, 
                event = trains$DiseaseSpecificDeath, 
                a=1, fold=10, oversample = T)
 
e2 <- glmnetboot(B=1000, X=chop, 
                 time = trains$DiseaseSpecificSurvival, 
                 event = trains$DiseaseSpecificDeath, 
                 a=0.5, fold=10, oversample = T)

e4 <- glmnetboot(B=1000, X=chop, 
                 time = trains$DiseaseSpecificSurvival, 
                 event = trains$DiseaseSpecificDeath, 
                 a=0.9, fold=10, oversample = T)

e5 <- glmnetboot(B=1000, X=chop, 
                 time = trains$DiseaseSpecificSurvival, 
                 event = trains$DiseaseSpecificDeath, 
                 a=0.9, fold=10, oversample = T)

r3 <- glmnetboot(B=1000, X=chop, 
                 time = trains$DiseaseSpecificSurvival, 
                 event = trains$DiseaseSpecificDeath, 
                 a=0, fold=10, oversample = F)

l3 <- glmnetboot(B=1000, X=chop, 
                 time = trains$DiseaseSpecificSurvival, 
                 event = trains$DiseaseSpecificDeath, 
                 a=1, fold=10, oversample = F)

e3 <- glmnetboot(B=1000, X=chop, 
                 time = trains$DiseaseSpecificSurvival, 
                 event = trains$DiseaseSpecificDeath, 
                 a=0.9, fold=10, oversample = F)

write.csv(r2$beta, file="ridgebeta2.csv")
write.csv(l2$beta, file="lassobeta2.csv")
write.csv(e2$beta, file="enetbeta2.csv")

write.csv(r3$beta, file="enetbeta3.csv")
write.csv(l3$beta, file="enetbeta3.csv")
write.csv(e3$beta, file="enetbeta3.csv")


# boxplot and barplot ===========
par(mar = c(7, 4.1, 2, 2.1))
boxplot(r2$beta, las=2,cex.axis=0.3)


boxplot(l2$beta, las=2,cex.axis=0.3)
barplot(l2$probzero, las=2, cex.names=0.3)

boxplot(e2$beta, las=2,cex.axis=0.3)
barplot(e2$probzero, las=2, cex.names=0.3)

# without oversample
boxplot(r3$beta, las=2,cex.axis=0.3)
barplot(r3$probzero, las=2, cex.names=0.3)

boxplot(l3$beta, las=2,cex.axis=0.3)
barplot(l3$probzero, las=2, cex.names=0.3)

boxplot(e3$beta, las=2,cex.axis=0.3)
barplot(e3$probzero, las=2, cex.names=0.3)

# ggplot
# ggbar
ggplot(data=as.data.frame(l2$probzero), aes(x= rownames(as.data.frame(l2$probzero)),
                                            y=l2$probzero)) + 
  geom_bar(stat = "identity") +theme(axis.text.x = element_text(size = rel(0.8), 
                                               angle = 90, hjust = 1)) + 
  ylab("probability of zero") + xlab("variable")

ggplot(data=as.data.frame(e4$probzero), aes(x= names(e4$probzero),y=e4$probzero)) + 
  geom_bar(stat = "identity") +theme(axis.text.x = element_text(size = rel(0.8), 
                                                                angle = 90, hjust = 1)) + 
  ylab("probability of zero") + xlab("variable")

ggplot(data=as.data.frame(l3$probzero), aes(x= rownames(as.data.frame(l2$probzero)),
                                            y=l2$probzero)) + 
  geom_bar(stat = "identity") +theme(axis.text.x = element_text(size = rel(0.8), 
                                                                angle = 90, hjust = 1)) + 
  ylab("probability of zero") + xlab("variable")

ggplot(data=as.data.frame(e3$probzero), aes(x= names(e4$probzero),y=e4$probzero)) + 
  geom_bar(stat = "identity") +theme(axis.text.x = element_text(size = rel(0.8), 
                                                                angle = 90, hjust = 1)) + 
  ylab("probability of zero") + xlab("variable")

# ggbox
ggplot(data=melt(as.data.frame(r2$beta)), aes(variable, value)) + 
  geom_boxplot() +theme(axis.text.x = element_text(size = rel(0.8), 
                                                   angle = 90, hjust = 1)) + 
  ylab("coefficient")
ggplot2::ggplot(r2$beta, aes(x= , y=))
labels<-names(chop)

ggplot(data=melt(as.data.frame(l2$beta)), aes(variable, value)) + 
  geom_boxplot() +theme(axis.text.x = element_text(size = rel(0.8), 
                                                   angle = 90, hjust = 1)) + 
  ylab("coefficient")

ggplot(data=melt(as.data.frame(e4$beta)), aes(variable, value)) + 
  geom_boxplot() +theme(axis.text.x = element_text(size = rel(0.8), 
                                                   angle = 90, hjust = 1)) + 
  ylab("coefficient")
# without oversample
ggplot(data=melt(as.data.frame(r3$beta)), aes(variable, value)) + 
  geom_boxplot() +theme(axis.text.x = element_text(size = rel(0.8), 
                                                   angle = 90, hjust = 1)) + 
  ylab("coefficient")

ggplot(data=melt(as.data.frame(l3$beta)), aes(variable, value)) + 
  geom_boxplot() +theme(axis.text.x = element_text(size = rel(0.8), 
                                                   angle = 90, hjust = 1)) + 
  ylab("coefficient")

ggplot(data=melt(as.data.frame(e3$beta)), aes(variable, value)) + 
  geom_boxplot() +theme(axis.text.x = element_text(size = rel(0.8), 
                                                   angle = 90, hjust = 1)) + 
  ylab("coefficient")
# CI and mean ==========
r2b <- r2$beta
r2b <- apply(r2b,2,sort,decreasing=F)
r2b <- data.frame(r2b)
# mean
barplot(abs(apply(r2b, 2, mean)), las=2, cex.names=0.3)
barplot(apply(r2b, 2, mean), las=2, cex.names=0.3)
# median 
barplot(abs((t(r2b))[,c(501)]), las=2, cex.names=0.3)
abline(h=.5)

r2t <- (t(r2b))[,c(25,975)]
r2t

rt2 <- t(r2t)
data.frame(rt2)$CD68.CD163.in.CT
rt2$
rt2 <- data.frame(t(r2t))
rt2 <- t(r2t)

# lasso
l2b <- l2$beta
l2b <- apply(l2b,2,sort,decreasing=F)
l2b <- data.frame(l2b)
# mean
barplot(abs(apply(l2b, 2, mean)), las=2, cex.names=0.3)
barplot(apply(l2b, 2, mean), las=2, cex.names=0.3)
# CI
l2t <- (t(l2b))[,c(25,975)]

l20 <- l2$probzero[order(l2$probzero)]
names(head(l20, 21))

# lasso oversampled
l3b <- l3$beta
l3b <- apply(l3b,2,sort,decreasing=F)
l3b <- data.frame(l3b)
# mean
barplot(abs(apply(l3b, 2, mean)), las=2, cex.names=0.3)
barplot(apply(l3b, 2, mean), las=2, cex.names=0.3)
a<-apply(l3b, 2, mean)
a[order(a)]
apply(l3b, 2, mean)

# CI
l3t <- (t(l3b))[,c(25,975)]
l3t
l30 <- l3$probzero[order(l3$probzero)]
l30
names(head(l30, 21))

# enet
e2b <- e2$beta
e2b <- apply(e2b,2,sort,decreasing=F)
e2b <- data.frame(e2b)
# mean
barplot(abs(apply(e2b, 2, mean)), las=2, cex.names=0.3)
barplot(apply(e2b, 2, mean), las=2, cex.names=0.3)

e2t <- (t(e2b))[,c(25,975)]

# alpha =0.9
boxplot(e3$beta, las=2,cex.axis=0.3)
barplot(e3$probzero, las=2, cex.names=0.3)

e3b <- e3$beta
e3b <- apply(e3b,2,sort,decreasing=F)
e3b <- data.frame(e3b)

e30 <- e3$probzero[order(e3$probzero)]
e30

# mean
barplot(abs(apply(e3b, 2, mean)), las=2, cex.names=0.3)
barplot(apply(e3b, 2, mean), las=2, cex.names=0.3)
apply(e3b, 2, mean)
e3t <- (t(e3b))[,c(25,975)]
e3t

# alpha =0.9 
boxplot(e4$beta, las=2,cex.axis=0.3)
barplot(e4$probzero, las=2, cex.names=0.3)
e4b <- e4$beta
e4b <- apply(e4b,2,sort,decreasing=F)
e4b <- data.frame(e4b)
# mean
barplot(abs(apply(e4b, 2, mean)), las=2, cex.names=0.3)
barplot(apply(e4b, 2, mean), las=2, cex.names=0.3)

e4t <- (t(e4b))[,c(25,975)]
e4t