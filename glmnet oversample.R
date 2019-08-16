# setwd("Data/")

# Output 
# overchop  - oversample design matrix with dummy variables for glmnet
# data4     - oversample data

# LOAD ======
# load(file = "Data/.RData")
# set seed ---------
set.seed(5099)

# time=brcancer[,1]
# status=brcancer[,2]
# geneexpr=as.matrix(brcancer[,-(1:2)])
# cox.lasso=cv.glmnet(geneexpr,Surv(time,status),family="cox",foldid=fold, standardize=FALSE)
# plot(cox.lasso)
# coefficients=coef(cox.lasso, s=cox.lasso$lambda.min)
# active.index=which(coefficients != 0)
# active.coefficients=coefficients[active.index]
# covarno=predict(cox.lasso, s=cox.lasso$lambda.min,type="nonzero")
# cbind(covarno,active.coefficients)

# You can control the randomness if you explicitly set foldid. Here an example for 5-fold CV
# library(caret)
# cvfold <- 5
# flds <- createFolds(data.sca$CaseNo, k = cvfold, list = TRUE, returnTrain = FALSE)
# foldids = rep(1,length(data.sca$CaseNo))
# foldids[flds$Fold2] = 2
# foldids[flds$Fold3] = 3
# foldids[flds$Fold4] = 4
# foldids[flds$Fold5] = 5
# Now run cv.glmnet with these foldids.
# lassoResults<-cv.glmnet(x=countDiffs,y=responseDiffs,alpha=1,foldid = foldids)
# load=========
save.image(file='oversample.RData')
load("oversample.RData") 

# parallel-----
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
# stopCluster(cl)

# library ====
library(ROSE)
library(DMwR)
library(gbm)
library(glmnet)
library(survival)
library(survminer)
# library(viridis)
# library(caret)
library(xtable)
library(stargazer)

# data input =====
source("data input.R")
source("rosesmoteover.R")

# oversampling ======
# train data with ovun.sample, ROSE and SMOTE
trains$DiseaseSpecificDeath <- as.factor(trains$DiseaseSpecificDeath)

data.ovun <- ovun.sample(DiseaseSpecificDeath ~ ., 
                                  data = trains, method = "both", 
                                  p=0.5, N = 1000, seed=5099)$data

table(data.ovun$DiseaseSpecificDeath)

data.rose <- ROSE(DiseaseSpecificDeath ~ ., data = trains, N=1000, seed = 5099)$data
table(data.rose$DiseaseSpecificDeath)

data.smote <- SMOTE(DiseaseSpecificDeath ~., data = trains, perc.over = 3200, perc.under = 105, k = 5)
as.data.frame(table(data.smote$DiseaseSpecificDeath))

# choose ovun.sample, rose or smote 
# data4 <- data.ovun/data.rose/data.smote
overchop <- data.smote[,-c(62:63)]
names(overchop)
# create dummy for glmnet
# xfactor <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=overchop)[, -1]
# # combine
# overchop <- as.matrix(data.frame(overchop[-c(62:67)], xfactor))
# overchop <- data.frame(overchop)
# chop <- data.sca[,-c(1,63:64)]

# for reproducible 10 fold cv.glmnet
# ind <- sample(1:113,113)
# data.sca <- data.sca[ind,]
# fold <- rep(1:10,length.out=113)

ind <- sample(1:nrow(overchop),nrow(overchop))
data.sca <- data.sca[ind,]
fold <- rep(1:10,length.out=nrow(overchop))


# change to dummy for glmnet
xfactor.o <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=overchop)[, -1]
overchop <- as.matrix(data.frame(overchop[-c(62:67)], xfactor.o))
overchop <- data.frame(overchop)

# xfactorv1 <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=chop.v1)[, -1]
# chop.v1 <- as.matrix(data.frame(chop.v1[-c(62:67)], xfactorv1))
# chop.v1 <- data.frame(chop.v1)
# 
# xfactorv2 <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=chop.v2)[,-1]
# chop.v2 <- as.matrix(data.frame(chop.v2[-c(62:67)], xfactorv2))
# chop.v2 <- data.frame(chop.v2)
# MODELLING ====
# https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#cox
# 5 fold or 10 fold
# ridge --------
cv.fit.r <- cv.glmnet(as.matrix(overchop), 
                    Surv(data4$DiseaseSpecificSurvival, data4$DiseaseSpecificDeath),
                    alpha = 0, family="cox", foldid=fold, maxit = 10000)
fit.r <- glmnet(as.matrix(overchop), 
              Surv(data4$DiseaseSpecificSurvival, data4$DiseaseSpecificDeath), 
              alpha = 0, family = "cox", maxit = 10000)
cv.fit.r$lambda.min
par(mfrow=c(1,2))
plot(fit.r, xvar="lambda")
abline(v=log(cv.fit.r$lambda.min))
title("ridge", line = -1)
plot(cv.fit.r)
abline(v=log(cv.fit.r$lambda.min))
title("ridge", line = -1)
cv.fit.r$lambda.min
pld.min.r <- cv.fit.r$cvm[cv.fit.r$lambda == cv.fit.r$lambda.min]

best_ridge_coef <- as.numeric(coef(cv.fit.r, s = cv.fit.r$lambda.min))

nr <- names(overchop[,Active.Index.r])

cr <- data.frame(nr, Active.Coefficients.r)
cr
cr <- cr[order(-abs(cr$Active.Coefficients.r)),]
cr
pr <- ggplot(cr, aes(x=seq_along(cr$nr), y=Active.Coefficients.r)) +
  geom_bar(stat="identity") + 
  labs(x= "index (ranked by coefficients)", y = "Coefficients of ridge regression")
pr
stargazer(cr)
cr


# best_ridge_coef
# LASSO ----------
cv.fit.l <- cv.glmnet(as.matrix(overchop), Surv(data4$DiseaseSpecificSurvival, 
                                        data4$DiseaseSpecificDeath), 
                    family="cox", foldid = fold, maxit = 2000000000)
fit.l <- glmnet(as.matrix(overchop), Surv(data4$DiseaseSpecificSurvival, 
                                  data4$DiseaseSpecificDeath), family = "cox", maxit = 2000000000)
plot_glmnet(fit.l, label=5)
par(mfrow=c(1,2))
plot_glmnet(fit.l)
vn=paste(names(overchop))
plot(fit.l)
vnat=coef(fit.l)
vnat=vnat[,ncol(vnat)] # remove the intercept, and get the coefficients at the end of the path
axis(4, at=vnat,line=-.5,label=vn,las=1,tick=FALSE, cex.axis=0.5) 

plot(fit.l, , xvar="lambda")
abline(v=log(cv.fit.l$lambda.min))
title("LASSO", line = -1)
plot(cv.fit.l)
abline(v=log(cv.fit.l$lambda.min))
title("LASSO", line = -1)
log(cv.fit.l$lambda.min)
cv.fit.l$lambda.min
pld.min.l <- cv.fit.l$cvm[cv.fit.l$lambda == cv.fit.l$lambda.min]
pld.min.l
Coefficients.l <- coef(fit.l, s = cv.fit.l$lambda.min)
Active.Index.l <- which(Coefficients.l != 0)
Active.Coefficients.l <- Coefficients.l[Active.Index.l]
Active.Coefficients.l
Active.Index.l
names(overchop[,Active.Index.l])
nl <- names(overchop[,Active.Index.l])

cl <- data.frame(nl, Active.Coefficients.l)
cl
cl <- cl[order(-abs(cl$Active.Coefficients.l)),]
cl
# elastic net ---------
# 0.5
cv.fit.e <- cv.glmnet(as.matrix(overchop), Surv(data4$DiseaseSpecificSurvival, 
                                            data4$DiseaseSpecificDeath), 
                      family="cox", foldid=fold, maxit = 2000000, alpha = 0.5)
fit.e <- glmnet(as.matrix(overchop), Surv(data4$DiseaseSpecificSurvival, 
                                      data4$DiseaseSpecificDeath), 
                family = "cox", alpha = 0.5)

plot(fit.e)
abline(v=log(cv.fit.e$lambda.min))
title("elastic net", line = -1)
plot(cv.fit.e)
abline(v=log(cv.fit.e$lambda.min))
title("elastic net", line = -1)
cv.fit.e$lambda.min

pld.min.e <- cv.fit.e$cvm[cv.fit.e$lambda == cv.fit.e$lambda.min]
pld.min.e

Coefficients.e <- coef(fit.e, s = cv.fit.e$lambda.min)
Active.Index.e <- which(Coefficients.e != 0)
Active.Coefficients.e <- Coefficients.e[Active.Index.e]
Active.Coefficients.e
Active.Index.e
names(overchop[,Active.Index.e])

# 0.9
cv.fit.e2 <- cv.glmnet(as.matrix(overchop), Surv(data4$DiseaseSpecificSurvival, 
                                                data4$DiseaseSpecificDeath), 
                      family="cox", foldid=fold, maxit = 2000000, alpha = 0.9)
fit.e2 <- glmnet(as.matrix(overchop), Surv(data4$DiseaseSpecificSurvival, 
                                          data4$DiseaseSpecificDeath), 
                family = "cox", alpha = 0.9)
fit.e2
plot(fit.e2, xvar="lambda")
abline(v=log(cv.fit.e2$lambda.min))
title("elastic net", line = -1)
plot(cv.fit.e2)
abline(v=log(cv.fit.e2$lambda.min))
title("elastic net", line = -1)
cv.fit.e2$lambda.min

pld.min.e2 <- cv.fit.e2$cvm[cv.fit.e2$lambda == cv.fit.e2$lambda.min]
pld.min.e2

Coefficients.e2 <- coef(fit.e2, s = cv.fit.e2$lambda.min)
Active.Index.e2 <- which(Coefficients.e2 != 0)
Active.Coefficients.e2 <- Coefficients.e2[Active.Index.e2]
Active.Coefficients.e2
Active.Index.e2
names(overchop[,Active.Index.e2])

ne2 <- names(overchop[,Active.Index.e2])

ce2 <- data.frame(ne2, Active.Coefficients.e2)
ce2
ce2 <- ce2[order(-abs(ce2$Active.Coefficients.e2)),]
ce2

# adaptive LASSO -------
# gamma = .5 #########
cv.fit.al.5 <- cv.glmnet(as.matrix(overchop),
                       Surv(data4$DiseaseSpecificSurvival, 
                            data4$DiseaseSpecificDeath), 
                       family="cox", foldid=fold, maxit=10000000, 
                       penalty.factor = 1 / abs(best_ridge_coef)^(.5))

fit.al.5 <- glmnet(as.matrix(overchop), Surv(data4$DiseaseSpecificSurvival, 
                                       data4$DiseaseSpecificDeath), 
                 family = "cox", penalty.factor = 1 / abs(best_ridge_coef)^(.5))
plot(fit.al.5)
abline(v=log(cv.fit.al.5$lambda.min))
title("Adaptive LASSO", line = -1)
plot(cv.fit.al.5)
abline(v=log(cv.fit.al.5$lambda.min))
title("Adaptive LASSO", line = -1)
cv.fit.al.5$lambda.min
pld.min.al.5 <- cv.fit.al.5$cvm[cv.fit.al.5$lambda == cv.fit.al.5$lambda.min]
pld.min.al.5

Coefficients.al.5 <- coef(fit.al.5, s = cv.fit.al.5$lambda.min)
Active.Index.al.5 <- which(Coefficients.al.5 != 0)
Active.Coefficients.al.5 <- Coefficients.al.5[Active.Index.al.5]
Active.Coefficients.al.5
Active.Index.al.5
names(overchop[,Active.Index.al.5])


# plot(cv.fit.al.5$glmnet.fit, xvar="lambda", label=TRUE)
# abline(v = log(cv.fit.al.5$lambda.min))
# abline(v = log(cv.fit.al.5$lambda.1se))
# coef(cv.fit.al.5, s=cv.fit.al.5$lambda.1se)
# coef <- coef(cv.fit.al.5, s='lambda.1se')
# selected_attributes <- (coef@i[-1]+1) ## Considering the structure of the data frame dataF as shown earlier

# gamma = 1 #########
cv.fit.al <- cv.glmnet(as.matrix(overchop),
                       Surv(data4$DiseaseSpecificSurvival, 
                                            data4$DiseaseSpecificDeath), 
                      family="cox", foldid=fold, maxit=10000000, 
                      penalty.factor = 1 / abs(best_ridge_coef), parallel = T)

fit.al <- glmnet(as.matrix(overchop), Surv(data4$DiseaseSpecificSurvival, 
                                      data4$DiseaseSpecificDeath), 
                family = "cox", penalty.factor = 1 / abs(best_ridge_coef))
plot(fit.al, xvar="lambda")
abline(v=log(cv.fit.al$lambda.min))
title("Adaptive LASSO", line = -1)
plot(cv.fit.al)
abline(v=log(cv.fit.al$lambda.min))
title("Adaptive LASSO", line = -1)
cv.fit.al$lambda.min
pld.min.al <- cv.fit.al$cvm[cv.fit.al$lambda == cv.fit.al$lambda.min]
pld.min.al
Coefficients.al <- coef(fit.al, s = cv.fit.al$lambda.min)
Active.Index.al <- which(Coefficients.al != 0)
Active.Coefficients.al <- Coefficients.al[Active.Index.al]
Active.Coefficients.al
Active.Index.al
names(overchop[,Active.Index.al])
nal <- names(overchop[,Active.Index.al])

cal <- data.frame(nal, Active.Coefficients.al)
cal
cal <- cal[order(-abs(cal$Active.Coefficients.al)),]
cal
# plotmo ==========
# plot(glmnet.cox)
# title("glmnet.cox", line=2)
# plot_glmnet(glmnet.cox, xvar="norm")
# plotres(glmnet.cox, which=3, do.par=FALSE)
# par(old.par)

# gamma = 2 to note ######## 
# bhat <- as.matrix(coef(cv.fit.l, s = "lambda.1se"))[-1, 1]
# if (all(bhat == 0)) {bhat <- rep(.Machine$double.eps * 2, length(bhat))
# }
# adpen <- (1/pmax(abs(bhat), .Machine$double.eps))
# m_adlasso <- glmnet(as.matrix(overchop),
#                     Surv(data4$DiseaseSpecificSurvival, 
#                          data4$DiseaseSpecificDeath), family = "cox", alpha = 1, 
#                     exclude = which(bhat ==0), penalty.factor = adpen)
# plot(m_adlasso)

cv.fit.al2 <- cv.glmnet(as.matrix(overchop),
                       Surv(data4$DiseaseSpecificSurvival,
                            data4$DiseaseSpecificDeath),
                       family="cox", foldid=fold,
                       penalty.factor = 1 / abs(best_ridge_coef)^2)

fit.al2 <- glmnet(as.matrix(overchop), Surv(data4$DiseaseSpecificSurvival,
                                       data4$DiseaseSpecificDeath),
                 family = "cox", penalty.factor = 1 / abs(best_ridge_coef)^2)


plot(fit.al2, xvar="lambda")
abline(v=log(cv.fit.al2$lambda.min))
title("Adaptive LASSO", line = -1)
plot(cv.fit.al2)
abline(v=log(cv.fit.al2$lambda.min))
title("Adaptive LASSO", line = -1)
cv.fit.al2$lambda.min
pld.min.al2 <- cv.fit.al2$cvm[cv.fit.al2$lambda == cv.fit.al2$lambda.min]
pld.min.al2
Coefficients.al2 <- coef(fit.al2, s = cv.fit.al2$lambda.min)
Active.Index.al2 <- which(Coefficients.al2 != 0)
Active.Coefficients.al2 <- Coefficients.al2[Active.Index.al2]
Active.Coefficients.al2
Active.Index.al2
names(overchop[,Active.Index.al2])

nal2 <- names(overchop[,Active.Index.al2])

cal2 <- data.frame(nal2, Active.Coefficients.al2)
cal2
cal2 <- cal2[order(-abs(cal2$Active.Coefficients.al2)),]
cal2

# gamma = 3 ##############
cv.fit.al3 <- cv.glmnet(as.matrix(overchop),
                        Surv(data4$DiseaseSpecificSurvival,
                             data4$DiseaseSpecificDeath),
                        family="cox", foldid=fold,
                        penalty.factor = 1 / abs(best_ridge_coef)^3)

fit.al3 <- glmnet(as.matrix(overchop), Surv(data4$DiseaseSpecificSurvival,
                                        data4$DiseaseSpecificDeath),
                  family = "cox", penalty.factor = 1 / abs(best_ridge_coef)^3)
plot(fit.al3, xvar="lambda")
abline(v=log(cv.fit.al3$lambda.min))
title("Adaptive LASSO", line = -1)
plot(cv.fit.al3)
abline(v=log(cv.fit.al3$lambda.min))
title("Adaptive LASSO", line = -1)
cv.fit.al3$lambda.min
pld.min.al3 <- cv.fit.al3$cvm[cv.fit.al3$lambda == cv.fit.al3$lambda.min]
pld.min.al3
Coefficients.al3 <- coef(fit.al3, s = cv.fit.al3$lambda.min)
Active.Index.al3 <- which(Coefficients.al3 != 0)
Active.Coefficients.al3 <- Coefficients.al3[Active.Index.al3]
Active.Coefficients.al3
Active.Index.al3
names(overchop[,Active.Index.al3])
nal3 <- names(overchop[,Active.Index.al3])

cal3 <- data.frame(nal3, Active.Coefficients.al3)
cal3
cal3 <- cal3[order(-abs(cal3$Active.Coefficients.al3)),]
cal3

# # ROC
# ## Extract predicted probabilities and observed outcomes.
# library(pROC)
# pY <- as.numeric(predict(fit.al, newx = chop, s = cv.fit.al$lambda.min, type = "response"))
# Y <- as.numeric()
# ## pROC for ROC construction
# roc1 <- pROC::roc(y ~ pY)
# ## Plot an ROC curve with AUC and threshold
# plot(roc1, print.auc = TRUE, print.thres = TRUE, print.thres.best.method = "youden")
# 
# SCAD ============
# cv.fit.s <- cv.ncvsurv(as.matrix(chop), Surv(data.sca$DiseaseSpecificSurvival, 
#                                             data.sca$DiseaseSpecificDeath), 
#                       penalty = "SCAD", seed = 5099, nfold = 10)
# scad_beta <- cv.fit.s$fit$beta[, cv.fit.s$min]
# model_compare[5, ] <- scad_beta[-1]

# # randomforest survival =============

data4$DiseaseSpecificDeath<-as.numeric(data4$DiseaseSpecificDeath)-1
data4$DiseaseSpecificDeath
v.obj <- rfsrc(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ .,
               data = data4,
               ntree = 5000, block.size = 10)
              


## print and plot the grow object
print(v.obj)
stargazer(v.obj)
plot(v.obj)
v.obj.imp <- vimp(v.obj)$importance
v.obj.imp.or <- v.obj.imp[order(abs(v.obj.imp))]
tail(v.obj.imp.or, 7)

plot(v.obj.imp)
abline(h=0.01)

## plot survival curves for first 10 individuals -- direct way
matplot(v.obj$time.interest, 100 * t(v.obj$survival.oob[1:10, ]),
        xlab = "Time", ylab = "Survival", type = "l", lty = 1)

## plot survival curves for first 10 individuals
## using function "plot.survival" 
plot.survival(v.obj, subset = 1:10)
p.v <- plot.variable(v.obj, surv.type = "surv", 
xvar.names = c("CD68.CD163..100.150.CD3", "CD68.100.150.CD3", 
               "CD68.150.200.CD8", "CD68..200.250.CD3", 
               "CD68.150.200.CD3", "CD68.CD163..150.200.CD3", 
               "CD68.CD163..200.250.CD3"),
                     partial = TRUE, smooth.lines = TRUE)
plot.variable(p.v)
p.v$plots.per.page <- 4
p.v$smooth.lines <- FALSE
plot.variable(p.v)

# prediction
# v1
rfv1 <- predict(v.obj, newdata=v1s)

v1s$rf <- ifelse(rfv1>0, "High" ,"Low") 

fitsurv.rf <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ rfv1, 
                        data=v1s, error="greenwood",conf.type="log") 
summary(fitsurv.rf)
ggsurvplot(fitsurv.rf,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ rfv1, 
         data=v1s)

fitsurv.rf <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ rfv1, 
                        data=v1s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ rfv1 ,
      data=v1s,ties="exact")

v1s$DiseaseSpecificDeath <- as.factor(v1s$DiseaseSpecificDeath)
v2s$DiseaseSpecificDeath <- as.factor(v2s$DiseaseSpecificDeath)

pred.test.fin <- predict(v.obj, 
                          newdata = v1s, 
                          importance = "none")
pred.test.fin2 <- predict(v.obj, 
                          newdata = v2s, 
                          importance = "none")
pred.test.fin$predicted
crf1 <- rcorr.cens(-pred.test.fin$predicted , 
                   Surv(v1s$DiseaseSpecificSurvival, v1s$DiseaseSpecificDeath))

crf2 <- rcorr.cens(-pred.test.fin2$predicted , 
                   Surv(v2s1$DiseaseSpecificSurvival, v2s1$DiseaseSpecificDeath))
crf <- cbind(crf1,crf2)

# gbm ==============
# data.sca1 <- na.omit(data.sca)
bfit <- gbm(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ .,
            data = data4, distribution = "coxph", 
            n.trees = 10000,
            interaction.depth = 1,
            shrinkage = 0.001,
            cv.folds = 5)
bfit2 <- gbm(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ .,
            data = data4, distribution = "coxph", 
            n.trees = 10000,
            interaction.depth = 2,
            shrinkage = 0.001,
            cv.folds = 5)
bfit2 <- gbm(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ .,
             data = data4, distribution = "coxph", 
             n.trees = 10000,
             interaction.depth = 3,
             shrinkage = 0.001,
             cv.folds = 5)
save.image(file='oversample.RData')
summary(bfit)
best.iter <- gbm.perf(bfit, method="cv")
plot(bfit)

plot(survfit(bfit), newdata=v1)

# create marginal plots
# plot variable X1,X2,X3 after "best" iterations
par(mfrow=c(1,3))
plot.gbm(bfit,1,best.iter)
plot.gbm(bfit,2,best.iter)
plot.gbm(bfit,3,best.iter)
par(mfrow=c(1,1))
plot.gbm(bfit,1:2,best.iter) # contour plot of variables 1 and 2 after "best" number iterations

# 3-way plots
plot.gbm(bfit,1:3,best.iter)

# print the first and last trees... just for curiosity
pretty.gbm.tree(bfit,1)
pretty.gbm.tree(bfit,bfit$n.trees)

# predict on the new data using "best" number of trees
# f.predict will be on the canonical scale (logit,log,etc.)
f.predict <- predict(gbm1,data2,best.iter)

# # Cox PH error
# # boosting
# risk <- rep(0,N)
# for(i in 1:N)
# {
#   risk[i] <- sum( (data2$tt>=data2$tt[i])*exp(f.predict) )
# }
# cat("Boosting:",sum( data2$delta*( f.predict - log(risk) ) ),"\n")
# 
# # linear model
# coxph1 <- coxph(Surv(tt,delta)~X1+X2+X3,data=data)
# f.predict <- predict(coxph1,newdata=data2)
# risk <- rep(0,N)
# for(i in 1:N)
# {
#   risk[i] <- sum( (data2$tt>=data2$tt[i])*exp(f.predict) )
# }
# cat("Linear model:",sum( data2$delta*( f.predict - log(risk) ) ),"\n")



# plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);
# plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
# file.copy(from=plots.png.paths, to="Plots/")

# prediction ========
# ridge
prv1 <- predict(object = cv.fit.r, newx = data.matrix(chop.v1), 
                s=cv.fit.r$lambda.min, type="response")

prv2 <- predict(object = cv.fit.r, newx = data.matrix(chop.v2), 
                s=cv.fit.r$lambda.min, type="response")

crv1 <- Hmisc::rcorr.cens(-prv1, Surv(v1s$DiseaseSpecificSurvival, 
                                      v1s$DiseaseSpecificDeath))

crv2 <- Hmisc::rcorr.cens(-prv2, Surv(v2s$DiseaseSpecificSurvival, 
                                      v2s$DiseaseSpecificDeath))
# lasso
plv1 <- predict(object = cv.fit.l, newx = data.matrix(chop.v1), 
                s=cv.fit.l$lambda.min, type="response")

plv2 <- predict(object = cv.fit.l, newx = data.matrix(chop.v2), 
                s=cv.fit.l$lambda.min, type="response")

clv1 <- Hmisc::rcorr.cens(-plv1, Surv(v1s$DiseaseSpecificSurvival, 
                                      v1s$DiseaseSpecificDeath))

clv2 <- Hmisc::rcorr.cens(-plv2, Surv(v2s$DiseaseSpecificSurvival, 
                                      v2s$DiseaseSpecificDeath))
# enet.5
pev1 <- predict(object = cv.fit.e, newx = data.matrix(chop.v1), 
                s=cv.fit.e$lambda.min, type="response")

pev2 <- predict(object = cv.fit.e, newx = data.matrix(chop.v2), 
                s=cv.fit.e$lambda.min, type="response")

cev1 <- Hmisc::rcorr.cens(-pev1, Surv(v1s$DiseaseSpecificSurvival, 
                                      v1s$DiseaseSpecificDeath))

cev2 <- Hmisc::rcorr.cens(-pev2, Surv(v2s$DiseaseSpecificSurvival, 
                                      v2s$DiseaseSpecificDeath))


# enet .9

pe2v1 <- predict(object = cv.fit.e2, newx = data.matrix(chop.v1), 
                 s=cv.fit.e2$lambda.min, type="response")

pe2v2 <- predict(object = cv.fit.e2, newx = data.matrix(chop.v2), 
                 s=cv.fit.e2$lambda.min, type="response")

ce2v1 <- Hmisc::rcorr.cens(-pe2v1, Surv(v1s$DiseaseSpecificSurvival, 
                                        v1s$DiseaseSpecificDeath))

ce2v2 <- Hmisc::rcorr.cens(-pe2v2, Surv(v2s$DiseaseSpecificSurvival, 
                                        v2s$DiseaseSpecificDeath))

# al5
pal5v1 <- predict(object = cv.fit.al.5, newx = data.matrix(chop.v1), 
                  s=cv.fit.al.5$lambda.min, type="response")

pal5v2 <- predict(object = cv.fit.al.5, newx = data.matrix(chop.v2), 
                  s=cv.fit.al.5$lambda.min, type="response")

cal5v1 <- Hmisc::rcorr.cens(-pal5v1, Surv(v1s$DiseaseSpecificSurvival, 
                                          v1s$DiseaseSpecificDeath))

cal5v2 <-Hmisc::rcorr.cens(-pal5v2, Surv(v2s$DiseaseSpecificSurvival, 
                                         v2s$DiseaseSpecificDeath))

# al
palv1 <- predict(object = cv.fit.al, newx = data.matrix(chop.v1), 
                 s=cv.fit.al$lambda.min, type="response")

palv2 <- predict(object = cv.fit.al, newx = data.matrix(chop.v2), 
                 s=cv.fit.al$lambda.min, type="response")

calv1 <- Hmisc::rcorr.cens(-palv1, Surv(v1s$DiseaseSpecificSurvival, 
                                        v1s$DiseaseSpecificDeath))

calv2 <- Hmisc::rcorr.cens(-palv2, Surv(v2s$DiseaseSpecificSurvival, 
                                        v2s$DiseaseSpecificDeath))


# al2
pal2v1 <- predict(object = cv.fit.al2, newx = data.matrix(chop.v1), 
                 s=cv.fit.al2$lambda.min, type="response")

pal2v2 <- predict(object = cv.fit.al2, newx = data.matrix(chop.v2), 
                 s=cv.fit.al2$lambda.min, type="response")

cal2v1 <- Hmisc::rcorr.cens(-pal2v1, Surv(v1s$DiseaseSpecificSurvival, 
                                        v1s$DiseaseSpecificDeath))

cal2v2 <- Hmisc::rcorr.cens(-pal2v2, Surv(v2s$DiseaseSpecificSurvival, 
                                        v2s$DiseaseSpecificDeath))

# al3
pal3v1 <- predict(object = cv.fit.al3, newx = data.matrix(chop.v1), 
                 s=cv.fit.al3$lambda.min, type="response")

pal3v2 <- predict(object = cv.fit.al3, newx = data.matrix(chop.v2), 
                 s=cv.fit.al3$lambda.min, type="response")

cal3v1 <- Hmisc::rcorr.cens(-pal3v1, Surv(v1s$DiseaseSpecificSurvival, 
                                        v1s$DiseaseSpecificDeath))

cal3v2 <- Hmisc::rcorr.cens(-pal3v2, Surv(v2s$DiseaseSpecificSurvival, 
                                        v2s$DiseaseSpecificDeath))

tt1 <- data.frame(crv1, clv1, ce2v1, calv1, cal2v1, cal3v1)
tt2 <- data.frame(crv2, clv2, ce2v2, calv2, cal2v2, cal3v2)
tt <- t(data.frame(tt1,tt2))
tt
stargazer::stargazer(tt)
# Surv =====
v1s$DiseaseSpecificDeath <- as.numeric(v1s$DiseaseSpecificDeath)
v2s$DiseaseSpecificDeath <- as.numeric(v2s$DiseaseSpecificDeath)
# ridge ------
# v1
r.predv1 <- predict(cv.fit.r, newx=as.matrix(chop.v1),s=c("lambda.min"), type="response")
v1s$r.min <- ifelse(r.predv1>1, "High" ,"Low")
v1s$r.min <- ifelse(scale(r.predv1 [,1], center = T, scale = T) > 0, "High" ,"Low") 
v1s$r.min
r.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ r.min , 
                        data=v1s, error="greenwood",conf.type="log") 
summary(r.fitsurv)
ggsurvplot(r.fitsurv,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ r.min, 
         data=v1s)

r.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ r.min, 
                        data=v1s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ r.min ,
      data=v1s,ties="exact")
# v2
r.predv2 <- predict(cv.fit.r, newx=as.matrix(chop.v2),s=c("lambda.min"), type="response")
v2s$r.min <- ifelse(r.predv2 [,1]>1, "High" ,"Low") 

r.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                          DiseaseSpecificDeath) ~ r.min , 
                     data=v2s, error="greenwood",conf.type="log") 
summary(r.fitsurv)
ggsurvplot(r.fitsurv,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ r.min, 
         data=v2s)

r.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                          DiseaseSpecificDeath) ~ r.min, 
                     data=v2s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ r.min ,
      data=v2s,ties="exact")
# lasso ========
l.predv1 <- predict(cv.fit.l, newx=as.matrix(chop.v1),s=c("lambda.min"), type="response")
v1s$l.min <- ifelse(l.predv1 [,1]>1, "High" ,"Low") 
l.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                          DiseaseSpecificDeath) ~ l.min , 
                     data=v1s, error="greenwood",conf.type="log") 
summary(l.fitsurv)
ggsurvplot(l.fitsurv,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ l.min, 
         data=v1s)

l.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                          DiseaseSpecificDeath) ~ l.min, 
                     data=v1s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ l.min ,
      data=v1s,ties="exact")
# v2
l.predv2 <- predict(cv.fit.l, newx=as.matrix(chop.v2),s=c("lambda.min"), type="response")
v2s$l.min <- ifelse(l.predv2 [,1]>10, "High" ,"Low") 
l.predv2
l.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                          DiseaseSpecificDeath) ~ l.min , 
                     data=v2s, error="greenwood",conf.type="log") 
summary(l.fitsurv)
ggsurvplot(l.fitsurv,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ l.min, 
         data=v2s)

l.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                          DiseaseSpecificDeath) ~ l.min, 
                     data=v2s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ l.min,
      data=v2s,ties="exact")
# adaptive 0.5
al.5.predv1 <- predict(cv.fit.al.5, newx=as.matrix(chop.v1),s=c("lambda.min"), type="response")
v1s$al.5.min <- ifelse(al.5.predv1 [,1]>1, "High" ,"Low") 
al.5.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                          DiseaseSpecificDeath) ~ al.5.min , 
                     data=v1s, error="greenwood",conf.type="log") 
summary(al.5.fitsurv)
ggsurvplot(al.5.fitsurv,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ al.5.min, 
         data=v1s)

al.5.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                          DiseaseSpecificDeath) ~ al.5.min, 
                     data=v1s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ al.5.min ,
      data=v1s,ties="exact")
# v2
al.5.predv2 <- predict(cv.fit.al.5, newx=as.matrix(chop.v2),s=c("lambda.min"), type="response")
v2s$al.5.min <- ifelse(al.5.predv2 [,1]>10, "High" ,"Low") 
al.5.predv2
al.5.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                          DiseaseSpecificDeath) ~ al.5.min , 
                     data=v2s, error="greenwood",conf.type="log") 
summary(al.5.fitsurv)
ggsurvplot(al.5.fitsurv,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ al.5.min, 
         data=v2s)

al.5.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                          DiseaseSpecificDeath) ~ al.5.min, 
                     data=v2s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ al.5.min,
      data=v2s,ties="exact")

# adaptive 1
al.predv1 <- predict(cv.fit.al, newx=as.matrix(chop.v1),s=c("lambda.min"), type="response")
v1s$al.min <- ifelse(al.predv1 [,1]>1, "High" ,"Low") 
al.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ al.min , 
                        data=v1s, error="greenwood",conf.type="log") 
summary(al.fitsurv)
ggsurvplot(al.fitsurv,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ al.min, 
         data=v1s)

al.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ al.min, 
                        data=v1s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ al.min ,
      data=v1s,ties="exact")
# v2
al.predv2 <- predict(cv.fit.al, newx=as.matrix(chop.v2),s=c("lambda.min"), type="response")
v2s$al.min <- ifelse(al.predv2 [,1]>1, "High" ,"Low") 
al.predv2
al.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ al.min , 
                        data=v2s, error="greenwood",conf.type="log") 
summary(al.fitsurv)
ggsurvplot(al.fitsurv,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ al.min, 
         data=v2s)

al.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ al.min, 
                        data=v2s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ al.min,
      data=v2s,ties="exact")
# adaptive 0.5
al.5.predv1 <- predict(cv.fit.al.5, newx=as.matrix(chop.v1),s=c("lambda.min"), type="response")
v1s$al.5.min <- ifelse(al.5.predv1 [,1]>1, "High" ,"Low") 
al.5.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ al.5.min , 
                        data=v1s, error="greenwood",conf.type="log") 
summary(al.5.fitsurv)
ggsurvplot(al.5.fitsurv,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ al.5.min, 
         data=v1s)

al.5.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ al.5.min, 
                        data=v1s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ al.5.min ,
      data=v1s,ties="exact")
# v2
al.5.predv2 <- predict(cv.fit.al.5, newx=as.matrix(chop.v2),s=c("lambda.min"), type="response")
v2s$al.5.min <- ifelse(al.5.predv2 [,1]>1, "High" ,"Low") 
al.5.predv2
al.5.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ al.5.min , 
                        data=v2s, error="greenwood",conf.type="log") 
summary(al.5.fitsurv)
ggsurvplot(al.5.fitsurv,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ al.5.min, 
         data=v2s)

al.5.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ al.5.min, 
                        data=v2s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ al.5.min,
      data=v2s,ties="exact")

# sessioninfo =========
sessionInfo()
