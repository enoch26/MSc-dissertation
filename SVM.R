source("rosesmoteover.R")
load("SVM.RData") 
svmFit.all <- train(DiseaseSpecificDeath ~ .,
                data = data4,
                method = "svmRadial", 
                trControl = cv_5_grid, 
                tuneLength = 8,
                metric = "ROC")
varImp(svmFit.all)
confusionMatrix(predict(svmFit.all), data4$DiseaseSpecificDeath)
ggplot(varImp(svmFit)) + theme(text = element_text(size=5)) 

# forward selection
svmFit <- train(DiseaseSpecificDeath ~ CD68.150.200.CD3+
                  CD68..200.250.CD3+
                  CD68pCD163n.150.200.CD3+
                  CD68.100.150.CD3+
                  CD3CT+
                  CD68pCD163n.200.250.CD3+
                  CD68pCD163n.100.150.CD3+
                  CD68.150.200.CD8+
                  CD68pCD163n.150.200.CD8,
                data = data4,
                method = "svmRadial", 
                trControl = cv_5_grid, 
                tuneLength = 8,
                metric = "ROC")

confusionMatrix(predict(svmFit), data4$DiseaseSpecificDeath)
confusionMatrix(predict(svmFit, newdata=v1s), v1s$DiseaseSpecificDeath)
confusionMatrix(predict(svmFit, newdata=v2s), v2s$DiseaseSpecificDeath)

svmFitttt <- train(DiseaseSpecificDeath ~ CD68.150.200.CD3+CD68..200.250.CD3+
                     CD68pCD163n.150.200.CD3+CD68.100.150.CD3+
                     CD3CT,
                data = data4,
                method = "svmRadial", 
                trControl = cv_5_grid, 
                tuneLength = 8,
                metric = "ROC")
confusionMatrix(predict(svmFitttt), data4$DiseaseSpecificDeath)
confusionMatrix(predict(svmFitttt, newdata=v1s), v1s$DiseaseSpecificDeath)
confusionMatrix(predict(svmFitttt, newdata=v2s), v2s$DiseaseSpecificDeath)
svmFit

plot(svmFit)
library(e1071)
svmksvm <- svm(DiseaseSpecificDeath ~ CD68.150.200.CD3+
                 CD68..200.250.CD3+
                 CD68pCD163n.150.200.CD3+
                 CD68.100.150.CD3+
                 CD3CT+
                 CD68pCD163n.200.250.CD3+
                 CD68pCD163n.100.150.CD3+
                 CD68.150.200.CD8+
                 CD68pCD163n.150.200.CD8,
               data = data4, kernel="radial", cost=32)

plot(svmksvm, data4, CD68.150.200.CD3~CD3CT)
plot(svmksvm, data4, CD68.150.200.CD3 ~ CD68..200.250.CD3)

svmksvm <- ksvm(DiseaseSpecificDeath ~ CD68.150.200.CD3+
                  CD68..200.250.CD3+
                  CD68pCD163n.150.200.CD3+
                  CD68.100.150.CD3+
                  CD3CT+
                  CD68pCD163n.200.250.CD3+
                  CD68pCD163n.100.150.CD3+
                  CD68.150.200.CD8+
                  CD68pCD163n.150.200.CD8,
                data = data4, type= "C-svc", kernel="rbfdot", C=32)

plot(svmksvm, )

# make a grid of the predictors
# rng11 <- range(data4$CD68.150.200.CD3)
# rng29 <- range(data4$CD68pCD163n.150.200.CD3)
# 
# n_points <- 100
# grid <- expand.grid(V11 = seq(rng11[1], rng11[2], length = n_points),
#                     V29 = seq(rng29[1], rng29[2], length = n_points))
# grid$decision <- predict(svmksvm, grid, type = "decision")[, 1]
# 
# 
# ggplot(data4, aes(x = V11, y = V29)) + 
#   geom_point(aes(col = Class)) +
#   geom_contour(data = grid, aes(z = decision), 
#                breaks = 0, col = "black") 


svmFitr1 <- train(DiseaseSpecificDeath ~ CD68pCD163n.overCD163.in.IM+
                    CD68.CD163.in.IM + Sex + Site+ pT + Diff + EMLVI,
                data = data4,
                method = "svmRadial", 
                trControl = cv_5_grid, 
                tuneLength = 8,
                metric = "ROC")
svmFitr1
varImp(svmFitr1)
confusionMatrix(predict(svmFitr1), data4$DiseaseSpecificDeath)
confusionMatrix(predict(svmFitr1, newdata=v1s), v1s$DiseaseSpecificDeath)
confusionMatrix(predict(svmFitr1, newdata=v2s), v2s$DiseaseSpecificDeath)

summary(svmFit)
svmvarimp <- varImp(svmFit)

svmvarimp

levels(v1s$DiseaseSpecificDeath) <- c("Censored", "Died")
make.names(levels(v1s$DiseaseSpecificDeath))
levels(v2s$DiseaseSpecificDeath) <- c("Censored", "Died")
make.names(levels(v2s$DiseaseSpecificDeath))

varImp(svmFit)
confusionMatrix(predict(svmFit), data4$DiseaseSpecificDeath)
confusionMatrix(predict(svmFit, newdata=v1s), v1s$DiseaseSpecificDeath)
confusionMatrix(predict(svmFit, newdata=v2s), v2s$DiseaseSpecificDeath)

# v1
svm.predv1 <- predict(svmFit, newdata=v1s)
v1s$svm<- ifelse(as.numeric(svm.predv1)>1, "High" ,"Low")
v1s$svm
svm.fitsurv1 <- survfit(Surv(DiseaseSpecificSurvival,
                             as.numeric(DiseaseSpecificDeath)) ~ svm , 
                        data=v1s, error="greenwood",conf.type="log") 
summary(svm.fitsurv1)
ggsurvplot(svm.fitsurv1,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ svm, 
         data=v1s)

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ svm,
      data=v1s,ties="exact")

# v2
svm.predv2 <- predict(svmFit, newdata=v2s)
v2s$svm<- ifelse(as.numeric(svm.predv2)>1, "High" ,"Low")
v2s$svm
svm.fitsurv2 <- survfit(Surv(DiseaseSpecificSurvival,
                             as.numeric(DiseaseSpecificDeath)) ~ svm , 
                        data=v2s, error="greenwood",conf.type="log") 
summary(svm.fitsurv2)
ggsurvplot(svm.fitsurv2,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ svm, 
         data=v2s)

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ svm,
      data=v2s,ties="exact")

# oversample
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

svmFitdat1 <- train(DiseaseSpecificDeath ~ CD68.150.200.CD3+
                      CD68..200.250.CD3+
                      CD68pCD163n.150.200.CD3+
                      CD68.100.150.CD3+
                      CD3CT+
                      CD68pCD163n.200.250.CD3+
                      CD68pCD163n.100.150.CD3+
                      CD68.150.200.CD8+
                      CD68pCD163n.150.200.CD8, 
                   data = data5Train[-62],
                   method = "svmRadial", 
                   trControl = cv_5_grid, 
                   tuneLength = 8,
                   metric = "ROC")

plot(svmFitdat1$finalModel, data5Train, DiseaseSpecificDeath ~ CD68.150.200.CD3)


# svmFitdat5 <- train(DiseaseSpecificDeath ~ CD68.150.200.CD3,
#                     data = v2,
#                     method = "svmRadial", 
#                     trControl = cv_5_grid, 
#                     tuneLength = 8,
#                     metric = "ROC")
# kernlab::plot(svmFitdat5$finalModel)


svmlight(DiseaseSpecificDeath~CD68.150.200.CD3+
           CD68..200.250.CD3+
           CD68pCD163n.150.200.CD3+
           CD68.100.150.CD3+
           CD3CT+
           CD68pCD163n.200.250.CD3+
           CD68pCD163n.100.150.CD3+
           CD68.150.200.CD8+
           CD68pCD163n.150.200.CD8, data=data5Train)

# partimat(DiseaseSpecificDeath~CD68.150.200.CD3+
#            CD68..200.250.CD3+
#            CD68pCD163n.150.200.CD3+
#            CD68.100.150.CD3+
#            CD3CT+
#            CD68pCD163n.200.250.CD3+
#            CD68pCD163n.100.150.CD3+
#            CD68.150.200.CD8+
#            CD68pCD163n.150.200.CD8, data=data5Train, method="svmlight",
#            plot.matrix = TRUE, imageplot = TRUE)

svmFitdat1
plot(svmFitdat1$finalModel, data = data5Test)
summary(svmFitdat1$finalModel)
class(svmFitdat1$finalModel)

confusionMatrix(predict(svmFitdat1), data5Train$DiseaseSpecificDeath)
confusionMatrix(predict(svmFitdat1, newdata=data5Test), data5Test$DiseaseSpecificDeath)

svm.predtest <- predict(svmFitdat1, newdata=data5Test)
data5Test$svm<- ifelse(as.numeric(svm.predtest)>1, "High" ,"Low")
data5Test$svm
as.numeric(data5Test$DiseaseSpecificDeath)

svm.fitsurtest <- survfit(Surv(DiseaseSpecificSurvival,
                              as.numeric(DiseaseSpecificDeath)) ~ svm , 
                         data=data5Test, error="greenwood",conf.type="log") 
summary(svm.fitsurtest)
ggsurvplot(svm.fitsurtest,data=data5Test,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ svm, 
         data=data5Test)

svm.cox <- coxph(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ svm,
      data=data5Test)

svm.ph <- cox.zph(svm.cox)

ggcoxzph(svm.ph)


save.image(file='SVM.RData')
