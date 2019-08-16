# pT ==========
a <- as.numeric(v1s$pT)
b<- as.numeric(v1s$DiseaseSpecificDeath)
confusionMatrix(as.factor(a),as.factor(b))
confusionMatrix(as.numeric(v1s$pT), as.numeric(v1s$DiseaseSpecificDeath))
confusionMatrix(as.factor(as.numeric(v2s$pT)), as.factor(as.numeric(v2s$DiseaseSpecificDeath)))
confusionMatrix(as.factor(as.numeric(data5Test$pT)), as.factor(as.numeric(data5Test$DiseaseSpecificDeath)))

pT.fitsurv1 <- survfit(Surv(DiseaseSpecificSurvival,
                            as.numeric(DiseaseSpecificDeath)) ~ pT , 
                       data=v1s, error="greenwood",conf.type="log") 
summary(pT.fitsurv1)
ggsurvplot(pT.fitsurv1,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ pT, 
         data=v1s)

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ pT,
      data=v1s,ties="exact")

pT.fitsurv2 <- survfit(Surv(DiseaseSpecificSurvival,
                            as.numeric(DiseaseSpecificDeath)) ~ pT , 
                       data=v2s, error="greenwood",conf.type="log") 
summary(pT.fitsurv2)
ggsurvplot(pT.fitsurv2,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ pT, 
         data=v2s)

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ pT,
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

pT.fitsurtest <- survfit(Surv(DiseaseSpecificSurvival,
                               as.numeric(DiseaseSpecificDeath)) ~ pT , 
                          data=data5Test, error="greenwood",conf.type="log") 
summary(pT.fitsurtest)
ggsurvplot(pT.fitsurtest,data=data5Test,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ pT, 
         data=data5Test)

# Diff =====
Diff.fitsurv1 <- survfit(Surv(DiseaseSpecificSurvival,
                            as.numeric(DiseaseSpecificDeath)) ~ Diff , 
                       data=v1s, error="greenwood",conf.type="log") 
summary(Diff.fitsurv1)
ggsurvplot(Diff.fitsurv1,data=v1s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ Diff, 
         data=v1s)

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ Diff,
      data=v1s,ties="exact")

Diff.fitsurv2 <- survfit(Surv(DiseaseSpecificSurvival,
                            as.numeric(DiseaseSpecificDeath)) ~ Diff , 
                       data=v2s, error="greenwood",conf.type="log") 
summary(Diff.fitsurv2)
ggsurvplot(Diff.fitsurv2,data=v2s,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ Diff, 
         data=v2s)

coxph(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ Diff,
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

Diff.fitsurtest <- survfit(Surv(DiseaseSpecificSurvival,
                              as.numeric(DiseaseSpecificDeath)) ~ Diff , 
                         data=data5Test, error="greenwood",conf.type="log") 
summary(Diff.fitsurtest)
ggsurvplot(Diff.fitsurtest,data=data5Test,risk.table = T,break.x.by=12,
           surv.median.line = "hv",pval=T, conf.int = T)

survdiff(Surv(DiseaseSpecificSurvival,as.numeric(DiseaseSpecificDeath)) ~ Diff, 
         data=data5Test)
