source("data input.R")
source("rosesmoteover.R")

save.image(file='gbm.RData')
load("gbm.RData") 

set.seed(5099)

library(gbm)
library(survival)
library(survminer)

# gbm ==============
# bfit <- gbm(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ .,
#             data = trains, distribution = "coxph", n.tree = 1000, cv.folds = 5,
#             bag.fraction = 1)
# 
# bfithappy <- gbm(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ .,
#             data = trains, distribution = "coxph", n.tree = 1000, cv.folds = 10,
#             bag.fraction = 1)
# 
# bfit5 <- gbm(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ .,
#             data = trains, distribution = "coxph", n.tree = 2000, cv.folds = 5,
#             bag.fraction = 1)
bfit <- gbm(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ .,
                        data = trains, distribution = "coxph", cv.folds=10)
            
bfit.int <- gbm(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ .,
            data = data4, distribution = "coxph", n.tree = 1000,
            shrinkage=0.1,
            interaction.depth=2,       # 1: additive model, 2: two-way interactions, etc
            # bag.fraction = 1,        # subsampling fraction, 0.5 is probably best
            # train.fraction = 0.5,
            cv.folds = 10,
            n.minobsinnode = 10,       # minimum total weight needed in each node
            keep.data = TRUE,
            verbose = FALSE) 

bfit.int3 <- gbm(Surv(DiseaseSpecificSurvival, 
                     as.numeric(DiseaseSpecificDeath)) ~ .,
                data = data4, distribution = "coxph", n.tree = 1000,
                shrinkage=0.1,
                interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                # bag.fraction = 1,        # subsampling fraction, 0.5 is probably best
                # train.fraction = 0.5,
                cv.folds = 5,
                n.minobsinnode = 10,       # minimum total weight needed in each node
                keep.data = TRUE,
                verbose = FALSE) 

bfit.add <- gbm(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ .,
              data = data4, distribution = "coxph", n.tree = 100,
              shrinkage=0.1,
              interaction.depth=1,       # 1: additive model, 2: two-way interactions, etc
              # bag.fraction = 1,        # subsampling fraction, 0.5 is probably best
              # train.fraction = 0.5,
              cv.folds = 10,
              n.minobsinnode = 10,       # minimum total weight needed in each node
              # keep.data = TRUE,
              verbose = FALSE) 

bfit.add2 <- gbm(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ .,
                data = data4, distribution = "coxph", n.tree = 1000,
                shrinkage=0.1,
                interaction.depth=1,       # 1: additive model, 2: two-way interactions, etc
                # bag.fraction = 1,        # subsampling fraction, 0.5 is probably best
                # train.fraction = 0.5,
                cv.folds = 10,
                n.minobsinnode = 10,       # minimum total weight needed in each node
                # keep.data = TRUE,
                verbose = FALSE)

summary(bfit.int3)

bfit.pred <- predict(bfit.int3, n.trees = 1000)
lambda0 <- basehaz.gbm(t = data4$DiseaseSpecificSurvival, 
                       delta = data4$DiseaseSpecificDeath,
                       t.eval =   sort(unique(data4$DiseaseSpecificSurvival)),
                       cumulative = FALSE, f.x = bfit.pred , smooth=TRUE)

lambda0 <- basehaz.gbm(t = v1s$DiseaseSpecificSurvival, 
                       delta = v1s$DiseaseSpecificDeath,
                       t.eval =   sort(unique(v1s$DiseaseSpecificSurvival)),
                       cumulative = FALSE, f.x = bfit.pred , smooth=TRUE)
lambda0

gbm.predv1 <- predict(bfit.int3, newdata = v1s, n.trees = 1000)

v1s$gbm <- ifelse(exp(gbm.predv1)>1, "High" ,"Low")
# v1s$gbm <- ifelse(scale(gbm.predv1 [,1], center = T, scale = T) > 0, "High" ,"Low") 
v1s$gbm
gbm.fitsurv1 <- survfit(Surv(DiseaseSpecificSurvival,
                          as.numeric(DiseaseSpecificDeath)) ~ gbm , 
                     data=v1s, error="greenwood",conf.type="log") 
summary(gbm.fitsurv1)
ggsurvplot(gbm.fitsurv1,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ gbm, 
         data=v1s)

# v2
gbm.predv2 <- predict(bfit.add2, newdata = v2s, n.trees = 1000)

v2s$gbm <- ifelse(exp(gbm.predv2)>1, "High" ,"Low")
 
v2s$gbm
gbm.fitsurv2 <- survfit(Surv(DiseaseSpecificSurvival,
                             as.numeric(DiseaseSpecificDeath)) ~ gbm , 
                        data=v2s, error="greenwood",conf.type="log") 
summary(gbm.fitsurv2)
ggsurvplot(gbm.fitsurv2,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ gbm, 
         data=v2s)






best.iter <- gbm.perf(bfit.add2 ,method="cv") 
print(best.iter)

bs <- summary(bfit)

xtable(bs)

# plot the performance
best.iter <- gbm.perf(bfit5,method="OOB")  # returns out-of-bag estimated best number of trees
print(best.iter)

best.iter <- gbm.perf(bfit5,method="cv") # returns test set estimate of best number of trees
print(best.iter)

best.iter <- gbm.perf(bfit5,method="test") # returns test set estimate of best number of trees
print(best.iter)

plot(bfit5)

plot(survfit(bfit), newdata=v1)

pred.train <- predict(bfit5, trains, n.trees = best.iter)
pred.test <- predict(bfit5, v1s, n.trees = best.iter)
pred.test2 <- predict(bfit5, v2s, n.trees = best.iter)
Hmisc::rcorr.cens(-pred.train, Surv(trains$DiseaseSpecificSurvival, 
                                    trains$DiseaseSpecificDeath))
cbv1 <- Hmisc::rcorr.cens(-pred.test, Surv(v1s$DiseaseSpecificSurvival, 
                                   v1s$DiseaseSpecificDeath))
cbv2 <- Hmisc::rcorr.cens(-pred.test2, Surv(v2s$DiseaseSpecificSurvival, 
                                    v2s$DiseaseSpecificDeath))
cb <- cbind(cbv1,cbv2)
stargazer(cb)

# plot variable influence

summary(bfit5,n.trees=1)         # based on the first tree

bb <- summary(bfit5,n.trees=best.iter) # based on the estimated best number of trees

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
pretty.gbm.tree(bfit5,1)
pretty.gbm.tree(bfit,bfit$n.trees)

# predict on the new data using "best" number of trees
# f.predict will be on the canonical scale (logit,log,etc.)
f.predict <- predict(bfit,v1,best.iter)


# Cox PH error
# boosting
risk <- rep(0,N)
for(i in 1:N)
{
  risk[i] <- sum( (data2$tt>=data2$tt[i])*exp(f.predict) )
}
cat("Boosting:",sum( data2$delta*( f.predict - log(risk) ) ),"\n")

# linear model
coxph1 <- coxph(Surv(tt,delta)~X1+X2+X3,data=data)
f.predict <- predict(coxph1,newdata=data2)
risk <- rep(0,N)
for(i in 1:N)
{
  risk[i] <- sum( (data2$tt>=data2$tt[i])*exp(f.predict) )
}
cat("Linear model:",sum( data2$delta*( f.predict - log(risk) ) ),"\n")


# save the model to disk
saveRDS(bfit5, "./bfit5.rds")

# load the model
super_model <- readRDS("./bfit5.rds")
