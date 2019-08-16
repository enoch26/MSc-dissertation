 # data input =====
source("rosesmoteover.R")
source("data input.R")
save.image(file='DAforrf.RData')
load("DAforrf.RData") 

library(doParallel)
cl <- makePSOCKcluster(3)
registerDoParallel(cl)

set.seed(5059)
# library ======
library(caret)
library(MASS)
library(dplyr)
library(mda)
library(penalizedLDA)
library(ggplot2)
library(kernlab)
library(klaR)
library(xtable)
library(stargazer)
library(biotools)
library(MVN)
library(survival)
library(survminer)



# https://github.com/topepo/caret/tree/master/RegressionTests/Code
# setting ===========
# grid search
cv_5_grid <- trainControl(method = "repeatedcv", number = 5, repeats = 10, 
                          classProbs = TRUE, summaryFunction = twoClassSummary,
                          allowParallel = TRUE, savePredictions = T)
# random search
cv_5_rand <- trainControl(method = "repeatedcv", number = 5, repeats = 10, 
                          classProbs = TRUE, summaryFunction = twoClassSummary,
                          allowParallel = TRUE, savePredictions = T,
                          search = "random")

# remove categorical var
# names(overchop)
# train1 <- overchop[, -c(116:126)]

# data3 without categorical variables
data3 <- data4[,-c(62,64:69)]

# MVN test
# mres <- mvn(data = data4, mvnTest = "mardia") 
# mres$multivariateNormality
# mres <- mvn(data = data4, mvnTest = "hz") 
# mres$multivariateNormality
# mres <- mvn(data = data4, mvnTest = "royston", 
#             univariateTest = "SW", desc = TRUE) 
# mres$multivariateNormality
# mres$univariateNormality
# mres <- mvn(data = data4, mvnTest = "dh") 
# mres$multivariateNormality
# mres <- mvn(data = data4, mvnTest = "energy") 
# mres$multivariateNormality

# fligner.test(DiseaseSpecificDeath ~ ., data = data4)
# logistic regression ===========
# log.fit <- glm(as.factor(DiseaseSpecificDeath) ~ .,
#                family = "binomial", data = data3, maxit = 1000)
# summary(log.fit)

# facotr 0 and 1 not accepted for discriminant analsis
levels(data3$DiseaseSpecificDeath) <- c("Censored", "Died")
make.names(levels(data3$DiseaseSpecificDeath))
levels(data4$DiseaseSpecificDeath) <- c("Censored", "Died")
make.names(levels(data4$DiseaseSpecificDeath))

trains$DiseaseSpecificDeath <- as.factor(trains$DiseaseSpecificDeath)
v1s$DiseaseSpecificDeath <- as.factor(v1s$DiseaseSpecificDeath)
v2s$DiseaseSpecificDeath <- as.factor(v2s$DiseaseSpecificDeath)
levels(trains$DiseaseSpecificDeath) <- c("Censored", "Died")
make.names(levels(trains$DiseaseSpecificDeath))
levels(v1s$DiseaseSpecificDeath) <- c("Censored", "Died")
make.names(levels(v1s$DiseaseSpecificDeath))
levels(v2s$DiseaseSpecificDeath) <- c("Censored", "Died")
make.names(levels(v2s$DiseaseSpecificDeath))
# RF =======

varImp(rfFit4$finalModel)
rfFit4$finalModel
summary(rfFit4$finalModel)

# oversample and standardize
rfFit.all <- train(DiseaseSpecificDeath ~ ., data = data3,
                 trControl = cv_5_grid, 
                 method = "rf", 
                 
                 ## Specify which metric to optimize
                 metric = "ROC",
                 verbose = FALSE)

# no oversample but standardize
rfFit.t <- train(DiseaseSpecificDeath ~ ., data = trains[,-62],
                 trControl = cv_5_grid, 
                 method = "rf", 
                 
                 ## Specify which metric to optimize
                 metric = "ROC",
                 verbose = FALSE)

# top 10 with oversampling 
rfFit5 <- train(DiseaseSpecificDeath ~ CD68pCD163n.150.200.CD8+
                   CD3CT+
                   CD68.150.200.CD3+
                   CD68pCD163n.150.200.CD3+
                   CD68.150.200.CD8
                 , data = data3,
                 trControl = cv_5_grid, 
                 method = "rf", 
                 metric = "ROC",
                 verbose = FALSE)

# top 10 with oversampling 
rfFit10 <- train(DiseaseSpecificDeath ~ CD68pCD163n.150.200.CD8+
                   CD3CT+
                   CD68.150.200.CD3+
                   CD68pCD163n.150.200.CD3+
                   CD68.150.200.CD8+
                   CD68..200.250.CD3+
                   CD3WTS+
                   CD3CD8.0.50.TB.TB.Number+
                   CD68pCD163n.200.250.CD8+
                   CD68pCD163n.100.150.CD3
                   , data = data3,
                  trControl = cv_5_grid, 
                  method = "rf", 
                  metric = "ROC",
                  verbose = FALSE)

rfFit20 <- train(DiseaseSpecificDeath ~ CD68pCD163n.150.200.CD8+
                   CD3CT+
                   CD68.150.200.CD3+
                   CD68pCD163n.150.200.CD3+
                   CD68.150.200.CD8+
                   CD68..200.250.CD3+
                   CD3WTS+
                   CD3CD8.0.50.TB.TB.Number+
                   CD68pCD163n.200.250.CD8+
                   CD68pCD163n.100.150.CD3+
                   CD68pCD163n.100.150.CD8+
                   CD68pCD163n.200.250.CD3+
                   CD68.100.150.CD8+
                   CD68.100.150.CD3+
                   CD68..200.250.CD8+
                   CD3IM+
                   CD163.200.250.CD3+
                   CD3.0.100.TB.TB.number+
                   CD68pCD163n.50.100.CD3+
                   CD8.0.50.TB.TB.number
                   , data = data3,
                  trControl = cv_5_grid, 
                  method = "rf", 
                  metric = "ROC",
                  verbose = FALSE)
# Varimp
summary(rfFit.all)
varImp(rfFit.all)
summary(rfFit.t)
summary(rfFit10)
plot(rfFit.all)
plot(rfFit.t)
plot(rfFit5)
plot(rfFit10)
plot(rfFit20)

# not bad with oversampling but with all
confusionMatrix(predict(rfFit.all), data3$DiseaseSpecificDeath)
confusionMatrix(predict(rfFit.all, newdata=v1s), v1s$DiseaseSpecificDeath)
confusionMatrix(predict(rfFit.all, newdata=v2s), v2s$DiseaseSpecificDeath)

# v1
rf.all.predv1 <- predict(rfFit.all, newdata=v1s)
v1s$rf<- ifelse(as.numeric(rf.all.predv1)>1, "High" ,"Low")
v1s$rf
rf.all.fitsurv1 <- survfit(Surv(DiseaseSpecificSurvival,
                            as.numeric(DiseaseSpecificDeath)) ~ rf , 
                       data=v1s, error="greenwood",conf.type="log") 
summary(rf.all.fitsurv1)
ggsurvplot(rf.all.fitsurv1,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ rf, 
         data=v1s)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ rf,
      data=v1s,ties="exact")

# v2
rf.all.predv2 <- predict(rfFit.all, newdata=v2s)
v2s$rf<- ifelse(as.numeric(rf.all.predv2)>1, "High" ,"Low")
v2s$rf
rf.all.fitsurv2 <- survfit(Surv(DiseaseSpecificSurvival,
                            as.numeric(DiseaseSpecificDeath)) ~ rf , 
                       data=v2s, error="greenwood",conf.type="log") 
summary(rf.all.fitsurv2)
ggsurvplot(rf.all.fitsurv2,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ rf, 
         data=v2s)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ rf,
      data=v2s,ties="exact")

# bad without oversampling
confusionMatrix(predict(rfFit.t), trains$DiseaseSpecificDeath)
confusionMatrix(predict(rfFit.t, newdata=v1s), v1s$DiseaseSpecificDeath)
confusionMatrix(predict(rfFit.t, newdata=v2s), v2s$DiseaseSpecificDeath)

# gets bad with only 5 predictors
confusionMatrix(predict(rfFit5), data3$DiseaseSpecificDeath)
confusionMatrix(predict(rfFit5, newdata=v1s), v1s$DiseaseSpecificDeath)
confusionMatrix(predict(rfFit5, newdata=v2s), v2s$DiseaseSpecificDeath)

# v1
gmb5.predv1 <- predict(rfFit5, newdata=v1s)
v1s$rf5<- ifelse(as.numeric(gmb5.predv1)>1, "High" ,"Low")
v1s$rf5
gmb5.fitsurv1 <- survfit(Surv(DiseaseSpecificSurvival,
                               as.numeric(DiseaseSpecificDeath)) ~ rf5 , 
                          data=v1s, error="greenwood",conf.type="log") 
summary(gmb5.fitsurv1)
ggsurvplot(gmb5.fitsurv1,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ rf5, 
         data=v1s)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ rf5,
      data=v1s,ties="exact")

# v2
gmb5.predv2 <- predict(rfFit5, newdata=v2s)
v2s$rf5<- ifelse(as.numeric(gmb5.predv2)>1, "High" ,"Low")
v2s$rf5
gmb5.fitsurv2 <- survfit(Surv(DiseaseSpecificSurvival,
                               as.numeric(DiseaseSpecificDeath)) ~ rf5 , 
                          data=v2s, error="greenwood",conf.type="log") 
summary(gmb5.fitsurv2)
ggsurvplot(gmb5.fitsurv2,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ rf5, 
         data=v2s)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ rf5,
      data=v2s,ties="exact")

# gets bad with only 10 predictors
confusionMatrix(predict(rfFit10), data3$DiseaseSpecificDeath)
confusionMatrix(predict(rfFit10, newdata=v1s), v1s$DiseaseSpecificDeath)
confusionMatrix(predict(rfFit10, newdata=v2s), v2s$DiseaseSpecificDeath)

# v1
gmb10.predv1 <- predict(rfFit10, newdata=v1s)
v1s$rf10<- ifelse(as.numeric(gmb10.predv1)>1, "High" ,"Low")
v1s$rf10
gmb10.fitsurv1 <- survfit(Surv(DiseaseSpecificSurvival,
                               as.numeric(DiseaseSpecificDeath)) ~ rf10 , 
                          data=v1s, error="greenwood",conf.type="log") 
summary(gmb10.fitsurv1)
ggsurvplot(gmb10.fitsurv1,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ rf10, 
         data=v1s)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ rf10,
      data=v1s,ties="exact")

# v2
gmb10.predv2 <- predict(rfFit10, newdata=v2s)
v2s$rf10<- ifelse(as.numeric(gmb10.predv2)>1, "High" ,"Low")
v2s$rf10
gmb10.fitsurv2 <- survfit(Surv(DiseaseSpecificSurvival,
                               as.numeric(DiseaseSpecificDeath)) ~ rf10 , 
                          data=v2s, error="greenwood",conf.type="log") 
summary(gmb10.fitsurv2)
ggsurvplot(gmb10.fitsurv2,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ rf10, 
         data=v2s)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ rf10,
      data=v2s,ties="exact")


# a bit better for 20
confusionMatrix(predict(rfFit20), data3$DiseaseSpecificDeath)
confusionMatrix(predict(rfFit20, newdata=v1s), v1s$DiseaseSpecificDeath)
confusionMatrix(predict(rfFit20, newdata=v2s), v2s$DiseaseSpecificDeath)

# v1
gmb20.predv1 <- predict(rfFit20, newdata=v1s)
v1s$rf20<- ifelse(as.numeric(gmb20.predv1)>1, "High" ,"Low")
v1s$rf20
gmb20.fitsurv1 <- survfit(Surv(DiseaseSpecificSurvival,
                                 as.numeric(DiseaseSpecificDeath)) ~ rf20 , 
                            data=v1s, error="greenwood",conf.type="log") 
summary(gmb20.fitsurv1)
ggsurvplot(gmb20.fitsurv1,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ rf20, 
         data=v1s)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ rf20,
      data=v1s,ties="exact")

# v2
gmb20.predv2 <- predict(rfFit20, newdata=v2s)
v2s$rf20<- ifelse(as.numeric(gmb20.predv2)>1, "High" ,"Low")
v2s$rf20
gmb20.fitsurv2 <- survfit(Surv(DiseaseSpecificSurvival,
                                 as.numeric(DiseaseSpecificDeath)) ~ rf20 , 
                            data=v2s, error="greenwood",conf.type="log") 
summary(gmb20.fitsurv2)
ggsurvplot(gmb20.fitsurv2,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ rf20, 
         data=v2s)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ rf20,
      data=v2s,ties="exact")



rfFit2$finalModel
summary(rfFit2$finalModel)

varImp(rfFit2$finalModel)
plot(rfFit2)

trellis.par.set(caretTheme())
plot(rfFit2, metric = "Kappa", plotType = "level",
     scales = list(x = list(rot = 90)))

ggplot(rfFit2) 

rfFit3 <- train(DiseaseSpecificDeath ~ ., data = data3,
                 trControl = cv_5_grid, 
                 method = "rf", 
                 
                 ## Specify which metric to optimize
                 metric = "ROC",
                 verbose = FALSE)
rfFit3

whichTwoPct <- tolerance(rfFit3$results, metric = "ROC", 
                         tol = 2, maximize = TRUE)  
cat("best model within 2 pct of best:\n")
rfFit3$results[whichTwoPct,1:7]

trellis.par.set(caretTheme())
plot(rfFit3) 
trellis.par.set(caretTheme())
plot(rfFit3, metric = "ROC", plotType = "level",
     scales = list(x = list(rot = 90)))
trellis.par.set(caretTheme())
densityplot(rfFit3, pch = "|")


# summary ==========
# https://topepo.github.io/caret/model-training-and-tuning.html
resamps <- resamples(list(rf = rfFit3,
                          SVM = svmFit,
                          RDA = fit_pda_grid))
resamps
summary(resamps)

theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
bwplot(resamps, layout = c(3, 1))

trellis.par.set(caretTheme())
dotplot(resamps, metric = "ROC")

trellis.par.set(theme1)
xyplot(resamps, what = "BlandAltman")

splom(resamps)

difValues <- diff(resamps)
difValues
summary(difValues)

trellis.par.set(theme1)
bwplot(difValues, layout = c(3, 1))

trellis.par.set(caretTheme())
dotplot(difValues)

# Surv

v1s$DiseaseSpecificDeath <- as.numeric(v1s$DiseaseSpecificDeath)
v2s$DiseaseSpecificDeath <- as.numeric(v2s$DiseaseSpecificDeath)

ldasurv1 <- survfit(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ lda, 
                    data=v2s, error = "greenwood", conf.type="log")

ldasurv1
ggsurvplot(ldasurv1,data=v1s,risk.table = T, surv.median.line = "hv", conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ lda, 
         data=v1s)

v1s$DiseaseSpecificDeath <- as.numeric(v1s$DiseaseSpecificDeath)
fitsurv <- survfit(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ svm, 
                   data=v1s, error="greenwood",conf.type="log")
summary(fitsurv)
ggsurvplot(fitsurv,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv")

svmsurv1 <- survfit(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ svm, 
                    data=v1s, error = "greenwood", conf.type = "log")

survdiff(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ svm, 
         data=v1s)


svmsurv2 <- survfit(Surv(DiseaseSpecificSurvival,DiseaseSpecificDeath) ~ svm, 
                    data=v2sna, conf.type = "log-log")
plot(svmsurv1)

summary(svmsurv1)
ggsurvplot(svmsurv1, data=v1s)
ggsurvplot(svmsurv, data=v2s, risk.table=T)

stopCluster(cl)