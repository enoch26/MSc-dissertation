source("rosesmoteover.R")
source("data input.R")

load("NB.RData") 
# Naive Bayes ============ 
# nbgrid <- data.frame(laplace=c(0,0.5,1.0), usekernel = TRUE, 
#                      adjust=c(0,0.5,1.0))
# 
# nbFit <- train(DiseaseSpecificDeath ~ ., data = data3,
#                method = "naive_bayes", 
#                trControl = cv_5_grid, 
#                tuneGrid = nbgrid,
#                metric = "ROC")

nbFit1 <- train(DiseaseSpecificDeath ~ .,
               data = data4,
               method = "naive_bayes", 
               trControl = cv_5_grid, 
               tuneGrid = nbgrid,
               metric = "ROC")
varImp(nbFit1)
confusionMatrix(predict(nbFit1), data4$DiseaseSpecificDeath)


nbFit <- train(DiseaseSpecificDeath ~ CD68.150.200.CD3+
                 CD68..200.250.CD3, 
               data = data4,
               method = "naive_bayes", 
               trControl = cv_5_grid, 
               tuneGrid = nbgrid,
               metric = "ROC")

nbFit
varImp(nbFit)
ggplot(varImp(nbFit)) + theme(text = element_text(size=5)) 

levels(v1s$DiseaseSpecificDeath) <- c("Censored", "Died")
make.names(levels(v1s$DiseaseSpecificDeath))
levels(v2s$DiseaseSpecificDeath) <- c("Censored", "Died")
make.names(levels(v2s$DiseaseSpecificDeath))

nbFit <- train(DiseaseSpecificDeath ~ CD68.150.200.CD3,
               data = data4,
               method = "naive_bayes", 
               trControl = cv_5_grid, 
               tuneGrid = nbgrid,
               metric = "ROC")

confusionMatrix(predict(nbFit), data4$DiseaseSpecificDeath)
confusionMatrix(predict(nbFit, newdata=v1s), v1s$DiseaseSpecificDeath)
confusionMatrix(predict(nbFit, newdata=v2s), v2s$DiseaseSpecificDeath)
# v1
nb.predv1 <- predict(nbFit, newdata=v1s)
v1s$nb<- ifelse(as.numeric(nb.predv1)>1, "High" ,"Low")
v1s$nb
nb.fitsurv1 <- survfit(Surv(DiseaseSpecificSurvival,
                            as.numeric(DiseaseSpecificDeath)) ~ nb , 
                       data=v1s, error="greenwood",conf.type="log") 
summary(nb.fitsurv1)
ggsurvplot(nb.fitsurv1,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ nb, 
         data=v1s)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ nb,
      data=v1s,ties="exact")

# v2
nb.predv2 <- predict(nbFit, newdata=v2s)
v2s$nb<- ifelse(as.numeric(nb.predv2)>1, "High" ,"Low")
v2s$nb
nb.fitsurv2 <- survfit(Surv(DiseaseSpecificSurvival,
                            as.numeric(DiseaseSpecificDeath)) ~ nb , 
                       data=v2s, error="greenwood",conf.type="log") 
summary(nb.fitsurv2)
ggsurvplot(nb.fitsurv2,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ nb, 
         data=v2s)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ nb,
      data=v2s,ties="exact")

# oversample
str(data1)

# data5 <- ROSE(DiseaseSpecificDeath ~ ., p=0.17, data = data1[,-c(1,126)], N=1000, seed = 5099)$data
# table(data5$DiseaseSpecificDeath)

# data5 <- ovun.sample(DiseaseSpecificDeath ~ .,
#                      data = data5, method = "both",
#                      p=0.17, N = 1000, seed=5099)$data

data5 <- SMOTE(DiseaseSpecificDeath ~., data = data1[,-c(1,126)], 
               perc.over = 390, perc.under = 710, k = 5)
as.data.frame(table(data5$DiseaseSpecificDeath))

trainIndex <- createDataPartition(data5$DiseaseSpecificDeath, p = .7, 
                                  list = FALSE, 
                                  times = 1)

levels(data5$DiseaseSpecificDeath) <- c("Censored", "Died")
make.names(levels(data5$DiseaseSpecificDeath))

data5Train <- data5[ trainIndex,]
data5Test  <- data5[-trainIndex,]

nbFitdat1 <- train(DiseaseSpecificDeath ~ CD68.150.200.CD3+
                 CD68..200.250.CD3, 
               data = data5Train[ ,-62],
               method = "naive_bayes", 
               trControl = cv_5_grid, 
               tuneGrid = nbgrid,
               metric = "ROC")

plot(nbFitdat1)

partimat(DiseaseSpecificDeath~CD68.150.200.CD3+
           CD68..200.250.CD3, data=data4, method="naiveBayes")

nbFitdat1$finalModel

confusionMatrix(predict(nbFitdat1), data5Train$DiseaseSpecificDeath)
confusionMatrix(predict(nbFitdat1, newdata=data5Test), data5Test$DiseaseSpecificDeath)

nb.predtest <- predict(nbFitdat1, newdata=data5Test)
data5Test$nb<- ifelse(as.numeric(nb.predtest)>1, "High" ,"Low")
data5Test$nb
as.numeric(data5Test$DiseaseSpecificDeath)

nb.fitsurtest <- survfit(Surv(DiseaseSpecificSurvival,
                            as.numeric(DiseaseSpecificDeath)) ~ nb , 
                       data=data5Test, error="greenwood",conf.type="log") 
summary(nb.fitsurtest)
ggsurvplot(nb.fitsurtest,data=data5Test,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ nb, 
         data=data5Test)

coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ nb,
      data=data5Test,ties="exact")

save.image(file='NB.RData')
