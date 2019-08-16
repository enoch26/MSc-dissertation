source("data input.R")


load("gbmno.RData") 

set.seed(5099)
library(dplyr)
library(ggplot2)
library(gbm)
library(survival)

# gbm ==============

bfit <- gbm(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ .,
            data = train, distribution = "coxph", n.tree = 100, cv.folds = 10)
# 
# bfithappy <- gbm(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ .,
#             data = trains, distribution = "coxph", n.tree = 1000, cv.folds = 10,
#             bag.fraction = 1)
# 
# bfit5 <- gbm(Surv(DiseaseSpecificSurvival, DiseaseSpecificDeath) ~ .,
#             data = trains, distribution = "coxph", n.tree = 2000, cv.folds = 5,
#             bag.fraction = 1)

bfit <- gbm(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ .,
            data = train, distribution = "coxph", n.tree = 10000,
            shrinkage=0.01,
            interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
            bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
            train.fraction = 0.5,
            cv.folds = 5,
            n.minobsinnode = 10,       # minimum total weight needed in each node
            keep.data = FALSE,
            verbose = FALSE)

bfit20 <- gbm(Surv(DiseaseSpecificSurvival, as.numeric(DiseaseSpecificDeath)) ~ .,
            data = train, distribution = "coxph", n.tree = 20000,
            shrinkage=0.1,
            interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
            bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
            train.fraction = 0.5,
            cv.folds = 5,
            n.minobsinnode = 10,       # minimum total weight needed in each node
            keep.data = FALSE,
            verbose = FALSE)


bfit

bfit2

bfit5

summary(bfit)
summary(bfit2)
summary(bfit5)

gbm.perf(bfit, 
         plot.it = TRUE, 
         oobag.curve = TRUE, 
         overlay = TRUE, 
         method = c("OOB","test")[1])

gbm.perf(bfit, 
         plot.it = TRUE, 
         oobag.curve = TRUE, 
         overlay = TRUE, 
         method = c("cv","test")[1])

gbm.perf(bfit2, 
         plot.it = TRUE, 
         oobag.curve = TRUE, 
         overlay = TRUE, 
         method = c("OOB","test")[1])

gbm.perf(bfit2, 
         plot.it = TRUE, 
         oobag.curve = TRUE, 
         overlay = TRUE, 
         method = c("cv","test")[1])

# select the best model
bfitt<-bfit

summary(bfitt)

Effects <- tibble::as_tibble(gbm::summary.gbm(bfitt, plotit = FALSE))
Effects %>% utils::head()
Effects %>% 
  # arrange descending to get the top influencers
  dplyr::arrange(desc(rel.inf)) %>%
  # sort to top 20
  dplyr::top_n(20) %>%
  # plot these data using columns
  ggplot(aes(x = forcats::fct_reorder(.f = var, 
                                      .x = rel.inf), 
             y = rel.inf, 
             fill = rel.inf)) +
  geom_col() +
  # flip
  coord_flip() +
  # format
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title = element_text()) + 
  xlab('variable') +
  ylab('relative influence')


gbm.predv1 <- predict(bfitt, newx=v1, type="response")

pred.train <- predict(bfitt, train, n.trees = best.iter)
time.interest <- sort(unique(train$DiseaseSpecificSurvival[train$DiseaseSpecificDeath==1]))
basehaz.cum  <- basehaz.gbm(t=train$DiseaseSpecificSurvival, 
                      delta=train$DiseaseSpecificDeath, 
                      pred.train,
                      t.eval=time.interest, 
                      cumulative = TRUE, 
                      smooth=T)
surf.i <- exp(-exp(pred.test[1])*basehaz.cum)
print(surf.i)

xtable(bs)

# plot the performance
best.iter <- gbm.perf(bfitt,method="OOB")  # returns out-of-bag estimated best number of trees
print(best.iter)

best.iter <- gbm.perf(bfitt,method="cv") # returns test set estimate of best number of trees
print(best.iter)

best.iter <- gbm.perf(bfitt,method="test") # returns test set estimate of best number of trees
print(best.iter)

par(mfrow=c(1,3))
a <- plot(bfitt, "CD68.100.150.CD3", best.iter, lwd = 2, col = "blue", ylab="measure of hazard ratio")
b <- plot(bfitt, "CD68.150.200.CD8", best.iter, lwd = 2, col = "blue", ylab="measure of hazard ratio")
c <- plot(bfitt, "CD3CT", best.iter, lwd = 2, col = "blue", ylab="measure of hazard ratio")
grid.arrange(a,b,c, ncol=2)


plot(bfitt, i.var = c("CD68.100.150.CD3","CD68.150.200.CD8","CD3CT"))

# i.var = 1, lwd = 2, col = "blue"

plot(survfit(bfit), newdata=v1)

pred.train <- predict(bfitt, trains, n.trees = best.iter)

pred.test <- predict(bfitt, v1s, n.trees = best.iter)
pred.test2 <- predict(bfitt, v2s, n.trees = best.iter)

# Hmisc::rcorr.cens(-pred.train, Surv(trains$DiseaseSpecificSurvival, 
#                                     trains$DiseaseSpecificDeath))
cbv1 <- Hmisc::rcorr.cens(-pred.test, Surv(v1$DiseaseSpecificSurvival, 
                                   as.numeric(v1$DiseaseSpecificDeath)))
cbv2 <- Hmisc::rcorr.cens(-pred.test2, Surv(v2$DiseaseSpecificSurvival, 
                                            as.numeric(v2$DiseaseSpecificDeath)))
cb <- cbind(cbv1,cbv2)
cb
stargazer(cb)

# plot variable influence

summary(bfitt,n.trees=1)         # based on the first tree

bb <- summary(bfitt,n.trees=best.iter) # based on the estimated best number of trees

# create marginal plots
# plot variable X1,X2,X3 after "best" iterations
par(mfrow=c(1,3))
plot.gbm(bfitt,1,best.iter)
plot.gbm(bfitt,2,best.iter)
plot.gbm(bfitt,3,best.iter)
par(mfrow=c(1,1))
plot.gbm(bfitt,1:2,best.iter) # contour plot of variables 1 and 2 after "best" number iterations

# 3-way plots
plot.gbm(bfitt,1:3,best.iter)

# print the first and last trees... just for curiosity
pretty.gbm.tree(bfitt,1)
pretty.gbm.tree(bfitt,bfitt$n.trees)

# predict on the new data using "best" number of trees
# f.predict will be on the canonical scale (logit,log,etc.)
f.predict <- predict(bfitt,v1,best.iter)


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

bfit5 <- super_model

save.image(file='gbmno.RData')
