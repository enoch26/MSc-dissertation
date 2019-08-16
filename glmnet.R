# setwd("Data/")

# set seed
set.seed(5099)

save.image(file="glmnet.RData")
load("glmnet.RData")

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


library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
# stopCluster(cl)

# library ====
library(dplyr)
library(survival)
library(glmnet)
# library(viridis)
# library(caret)
# library(xtable)
# library(stargazer)

# Data input ====
source("data input.R")

# change to dummy for glmnet
xfactor <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=chop)[, -1]
chop <- as.matrix(data.frame(chop[-c(62:67)], xfactor))
chop <- data.frame(chop)

xfactorv1 <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=chop.v1)[, -1]
chop.v1 <- as.matrix(data.frame(chop.v1[-c(62:67)], xfactorv1))
chop.v1 <- data.frame(chop.v1)

xfactorv2 <- model.matrix(~ Sex + Age + pT + Site + Diff + EMLVI , data=chop.v2)[,-1]
chop.v2 <- as.matrix(data.frame(chop.v2[-c(62:67)], xfactorv2))
chop.v2 <- data.frame(chop.v2)

# for reproducible 10 fold cv.glmnet
ind <- sample(1:113,113)
trains <- trains[ind,]
fold <- rep(1:10,length.out=113)

# MODELLING ====
# http://r.789695.n4.nabble.com/Interperting-results-of-glmnet-and-coxph-plot-Brier-score-and-Harrel-s-C-Index-am-I-doing-something--td4677166.html
# https://github.com/cran/FRESA.CAD/blob/master/R/crossValidationFeatureSelection.Res.R
# https://stats.stackexchange.com/questions/88696/discrepancy-between-log-likelihood-harrells-c-index-and-brier-score-for-the-eva
# https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#cox
# ridge --------
cv.fit.r <- cv.glmnet(data.matrix(chop), 
                    Surv(trains$DiseaseSpecificSurvival, trains$DiseaseSpecificDeath),
                    alpha = 0, family="cox", foldid=fold)
fit.r <- glmnet(as.matrix(chop), 
              Surv(trains$DiseaseSpecificSurvival, trains$DiseaseSpecificDeath), 
              alpha = 0, family = "cox")

fit.r$nobs
par(mfrow=c(1,2))
plot(fit.r)
abline(v=log(cv.fit.r$lambda.min))
log(cv.fit.r$lambda.min)
title("ridge", line = -1)
plot(cv.fit.r)
abline(v=log(cv.fit.r$lambda.min))
title("ridge", line = -1)
# optimal lambda
l.r <- cv.fit.r$lambda.min
l.r

# partial likelihood deviance
pld.min.r <- cv.fit.r$cvm[cv.fit.r$lambda == cv.fit.r$lambda.min]
pld.min.r

best_ridge_coef <- as.numeric(coef(cv.fit.r, s = cv.fit.r$lambda.min))
best_ridge_coef 

Coefficients.r <- coef(fit.r, s = cv.fit.r$lambda.min)
Active.Index.r <- which(Coefficients.r != 0)
Active.Coefficients.r <- Coefficients.r[Active.Index.r]
Active.Coefficients.r
Active.Index.r
nr <- names(chop[,Active.Index.r])
nr
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

# LASSO ----------
cv.fit.l <- cv.glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival, 
                                        trains$DiseaseSpecificDeath), 
                    family="cox", foldid = fold)
fit.l <- glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival, 
                                  trains$DiseaseSpecificDeath), family = "cox")
plot(fit.l)
abline(v=log(cv.fit.l$lambda.min))
title("LASSO", line = -1)
plot(cv.fit.l)
abline(v=log(cv.fit.l$lambda.min))
title("LASSO", line = -1)
log(cv.fit.l$lambda.min)
pld.min.l <- cv.fit.l$cvm[cv.fit.l$lambda == cv.fit.l$lambda.min]
pld.min.l
Coefficients.l <- coef(fit.l, s = cv.fit.l$lambda.min)
Active.Index.l <- which(Coefficients.l != 0)
Active.Coefficients.l <- Coefficients.l[Active.Index.l]
Active.Coefficients.l
Active.Index.l
nl <- names(chop[,Active.Index.l])
nl
cl <- data.frame(nl, Active.Coefficients.l)
stargazer(cl)
cl
cl <- cl[order(-abs(cl$Active.Coefficients.l)),]
cl
# elastic net ---------
# alpha =0.5 ##########
cv.fit.e <- cv.glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival, 
                                            trains$DiseaseSpecificDeath), 
                      family="cox", foldid=fold, alpha = 0.5)
fit.e <- glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival, 
                                      trains$DiseaseSpecificDeath), 
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
ne <- names(chop[,Active.Index.e])

ce <- data.frame(ne, Active.Coefficients.e)
ce
stargazer(ce)



# alpha = 0.9 #########
cv.fit.e2 <- cv.glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival, 
                                            trains$DiseaseSpecificDeath), 
                      family="cox", foldid=fold, alpha = 0.9)
fit.e2 <- glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival, 
                                      trains$DiseaseSpecificDeath), 
                family = "cox", alpha = 0.9)

plot(fit.e2)
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
ne2 <- names(chop[,Active.Index.e2])

ce2 <- data.frame(ne2, Active.Coefficients.e2)
ce2 <- ce2[order(-abs(ce2$Active.Coefficients.e2)),]
stargazer(ce2)
ce2


# adaptive LASSO -------
# gamma = .1 #########
cv.fit.al.1 <- cv.glmnet(as.matrix(chop),
                         Surv(trains$DiseaseSpecificSurvival, 
                              trains$DiseaseSpecificDeath), 
                         family="cox", foldid=fold, maxit=10000000, 
                         penalty.factor = 1 / abs(best_ridge_coef)^(.1))

fit.al.1 <- glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival, 
                                         trains$DiseaseSpecificDeath), 
                   family = "cox", penalty.factor = 1 / abs(best_ridge_coef))
plot(fit.al.1)
abline(v=log(cv.fit.al.1$lambda.min))
title("Adaptive LASSO", line = -1)
plot(cv.fit.al.1)
abline(v=log(cv.fit.al.1$lambda.min))
title("Adaptive LASSO", line = -1)
cv.fit.al.1$lambda.min
pld.min.al.1 <- cv.fit.al.1$cvm[cv.fit.al.1$lambda == cv.fit.al.1$lambda.min]
pld.min.al.1

Coefficients.al.1 <- coef(fit.al.1, s = cv.fit.al.1$lambda.min)
Active.Index.al.1 <- which(Coefficients.al.1 != 0)
Active.Coefficients.al.1 <- Coefficients.al.1[Active.Index.al.1]
Active.Coefficients.al.1
Active.Index.al.1
nal.1 <- names(chop[,Active.Index.al.1])

cal.1 <- data.frame(nal.1, Active.Coefficients.al.1)
stargazer(cal.1)
cal.1

pal.1v1 <- predict(object = cv.fit.al.1, newx = data.matrix(chop.v1), 
                  s=cv.fit.al.1$lambda.min, type="response")

pal.1v2 <- predict(object = cv.fit.al.1, newx = data.matrix(chop.v2), 
                  s=cv.fit.al.1$lambda.min, type="response")

cal.1v1 <- Hmisc::rcorr.cens(-pal.1v1, Surv(v1s$DiseaseSpecificSurvival, 
                                          v1s$DiseaseSpecificDeath))

cal.1v2 <-Hmisc::rcorr.cens(-pal.1v2, Surv(v2s$DiseaseSpecificSurvival, 
                                         v2s$DiseaseSpecificDeath))

cal.1v1
cal.1v2
# plot(cv.fit.al.5$glmnet.fit, xvar="lambda", label=TRUE)
# abline(v = log(cv.fit.al.5$lambda.min))
# abline(v = log(cv.fit.al.5$lambda.1se))
# coef(cv.fit.al.5, s=cv.fit.al.5$lambda.1se)
# coef <- coef(cv.fit.al.5, s='lambda.1se')
# selected_attributes <- (coef@i[-1]+1) ## Considering the structure of the trains frame trainsF as shown earlier

# gamma = .5 #########
cv.fit.al.5 <- cv.glmnet(as.matrix(chop),
                       Surv(trains$DiseaseSpecificSurvival, 
                            trains$DiseaseSpecificDeath), 
                       family="cox", foldid=fold, maxit=10000000, 
                       penalty.factor = 1 / abs(best_ridge_coef)^(.5))

fit.al.5 <- glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival, 
                                       trains$DiseaseSpecificDeath), 
                 family = "cox", penalty.factor = 1 / abs(best_ridge_coef))
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
nal5 <- names(chop[,Active.Index.al.5])

cal5 <- data.frame(nal5, Active.Coefficients.al.5)
cal5
cal5 <- cal5[order(-abs(cal5$Active.Coefficients.al.5)),]
stargazer(cal5)
cal5



# plot(cv.fit.al.5$glmnet.fit, xvar="lambda", label=TRUE)
# abline(v = log(cv.fit.al.5$lambda.min))
# abline(v = log(cv.fit.al.5$lambda.1se))
# coef(cv.fit.al.5, s=cv.fit.al.5$lambda.1se)
# coef <- coef(cv.fit.al.5, s='lambda.1se')
# selected_attributes <- (coef@i[-1]+1) ## Considering the structure of the trains frame trainsF as shown earlier

# gamma = 1 #########
cv.fit.al <- cv.glmnet(as.matrix(chop),
                       Surv(trains$DiseaseSpecificSurvival, 
                                            trains$DiseaseSpecificDeath), 
                      family="cox", foldid=fold, maxit=10000000, 
                      penalty.factor = 1 / abs(best_ridge_coef))

fit.al <- glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival, 
                                      trains$DiseaseSpecificDeath), 
                family = "cox", penalty.factor = 1 / abs(best_ridge_coef))
plot(fit.al)
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
nal <- names(chop[,Active.Index.al])

cal <- data.frame(nal, Active.Coefficients.al)
stargazer(cal)
cal

# plotmo ==========
# plot(glmnet.cox)
# title("glmnet.cox", line=2)
# plot_glmnet(glmnet.cox, xvar="norm")
# plotres(glmnet.cox, which=3, do.par=FALSE)
# par(old.par)

# gamma = 2 ######## 
# bhat <- as.matrix(coef(cv.fit.l, s = "lambda.1se"))[-1, 1]
# if (all(bhat == 0)) {bhat <- rep(.Machine$double.eps * 2, length(bhat))
# }
# adpen <- (1/pmax(abs(bhat), .Machine$double.eps))
# m_adlasso <- glmnet(as.matrix(chop),
#                     Surv(trains$DiseaseSpecificSurvival, 
#                          trains$DiseaseSpecificDeath), family = "cox", alpha = 1, 
#                     exclude = which(bhat ==0), penalty.factor = adpen)
# plot(m_adlasso)
# 
# cv.fit.al2 <- cv.glmnet(as.matrix(chop),
#                        Surv(trains$DiseaseSpecificSurvival,
#                             trains$DiseaseSpecificDeath),
#                        family="cox", foldid=fold,
#                        penalty.factor = 1 / abs(best_ridge_coef)^2)
# 
# fit.al2 <- glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival,
#                                        trains$DiseaseSpecificDeath),
#                  family = "cox", penalty.factor = 1 / abs(best_ridge_coef)^2)
# plot(fit.al2)
# abline(v=log(cv.fit.al2$lambda.min))
# title("Adaptive LASSO", line = -1)
# plot(cv.fit.al2)
# abline(v=log(cv.fit.al2$lambda.min))
# title("Adaptive LASSO", line = -1)
# cv.fit.al2$lambda.min
# pld.min.al2 <- cv.fit.al2$cvm[cv.fit.al2$lambda == cv.fit.al2$lambda.min]
# pld.min.al2
# Coefficients.al2 <- coef(fit.al2, s = cv.fit.al2$lambda.min)
# Active.Index.al2 <- which(Coefficients.al2 != 0)
# Active.Coefficients.al2 <- Coefficients.al2[Active.Index.al2]
# Active.Coefficients.al2
# Active.Index.al2
# nal2 <- names(chop[,Active.Index.al2])
# nal2
# # gamma = 3 ##############
# cv.fit.al3 <- cv.glmnet(as.matrix(chop),
#                         Surv(trains$DiseaseSpecificSurvival,
#                              trains$DiseaseSpecificDeath),
#                         family="cox", foldid=fold,
#                         penalty.factor = 1 / abs(best_ridge_coef)^3)
# 
# fit.al3 <- glmnet(as.matrix(chop), Surv(trains$DiseaseSpecificSurvival,
#                                         trains$DiseaseSpecificDeath),
#                   family = "cox", penalty.factor = 1 / abs(best_ridge_coef)^3)
# plot(fit.al3)
# abline(v=log(cv.fit.al3$lambda.min))
# title("Adaptive LASSO", line = -1)
# plot(cv.fit.al3)
# abline(v=log(cv.fit.al3$lambda.min))
# title("Adaptive LASSO", line = -1)
# cv.fit.al3$lambda.min
# pld.min.al3 <- cv.fit.al3$cvm[cv.fit.al3$lambda == cv.fit.al3$lambda.min]
# pld.min.al3
# Coefficients.al3 <- coef(fit.al3, s = cv.fit.al3$lambda.min)
# Active.Index.al3 <- which(Coefficients.al3 != 0)
# Active.Coefficients.al3 <- Coefficients.al3[Active.Index.al3]
# Active.Coefficients.al3
# Active.Index.al3
# names(chop[,Active.Index.al3])

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
# cv.fit.s <- cv.ncvsurv(as.matrix(chop), Surv(trains.sca$DiseaseSpecificSurvival, 
#                                             trains.sca$DiseaseSpecificDeath), 
#                       penalty = "SCAD", seed = 5099, nfold = 10)
# scad_beta <- cv.fit.s$fit$beta[, cv.fit.s$min]
# model_compare[5, ] <- scad_beta[-1]



# plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);
# plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
# file.copy(from=plots.png.paths, to="Plots/")

# Prediction =======
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
# pev1 <- predict(object = cv.fit.e, newx = data.matrix(chop.v1), 
#                 s=cv.fit.e$lambda.min, type="response")
# 
# pev2 <- predict(object = cv.fit.e, newx = data.matrix(chop.v2), 
#                 s=cv.fit.e$lambda.min, type="response")
# 
# cev1 <- Hmisc::rcorr.cens(-pev1, Surv(v1s$DiseaseSpecificSurvival, 
#                                       v1s$DiseaseSpecificDeath))
# 
# cev2 <- Hmisc::rcorr.cens(-pev2, Surv(v2s$DiseaseSpecificSurvival, 
#                                       v2s$DiseaseSpecificDeath))


# enet .9

pe2v1 <- predict(object = cv.fit.e2, newx = data.matrix(chop.v1), 
                 s=cv.fit.e2$lambda.min, type="response")

pe2v2 <- predict(object = cv.fit.e2, newx = data.matrix(chop.v2), 
                 s=cv.fit.e2$lambda.min, type="response")

ce2v1 <- Hmisc::rcorr.cens(-pe2v1, Surv(v1s$DiseaseSpecificSurvival, 
                                        v1s$DiseaseSpecificDeath))

ce2v2 <- Hmisc::rcorr.cens(-pe2v2, Surv(v2s$DiseaseSpecificSurvival, 
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

tt1 <- data.frame(crv1, clv1, cev1, calv1)
tt2 <- data.frame(crv2, clv2, cev2, calv2)
tt <- data.frame(tt1,tt2)
tt
stargazer::stargazer(tt)
# Surv ========
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
# lasso 
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
v2s$l.min <- ifelse(l.predv2 [,1]>1, "High" ,"Low") 
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
# elastic net alpha = 0.9
e2.predv1 <- predict(cv.fit.e2, newx=as.matrix(chop.v1),s=c("lambda.min"), type="response")
v1s$e2.min <- ifelse(e2.predv1 [,1]>1, "High" ,"Low") 
e2.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ e2.min , 
                        data=v1s, error="greenwood",conf.type="log") 
summary(e2.fitsurv)
ggsurvplot(e2.fitsurv,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ e2.min, 
         data=v1s)

e2.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ e2.min, 
                        data=v1s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ e2.min ,
      data=v1s,ties="exact")
# v2
e2.predv2 <- predict(cv.fit.e2, newx=as.matrix(chop.v2),s=c("lambda.min"), type="response")
v2s$e2.min <- ifelse(e2.predv2 [,1]>1, "High" ,"Low") 
e2.predv2
e2.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ e2.min , 
                        data=v2s, error="greenwood",conf.type="log") 
summary(e2.fitsurv)
ggsurvplot(e2.fitsurv,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ e2.min, 
         data=v2s)

e2.fitsurv <- survfit(Surv(DiseaseSpecificSurvival,
                             DiseaseSpecificDeath) ~ e2.min, 
                        data=v2s, error="greenwood",conf.type="log")

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ e2.min,
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

stargazer(tt)
# sessioninfo =========
sessionInfo()
