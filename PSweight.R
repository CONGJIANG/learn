library(PSweight)

data("psdata")
ps.formula <- trt ~ cov1 + cov2+cov3+cov4+cov5+cov6
msstat <- SumStat(ps.formula, trtgrp="2", data=psdata,
                  weight=c("ATE","ATO","ATT"))
plot(msstat, type="hist")
plot(msstat, type="balance", weighted.var=TRUE, threshold=0.1, metric="ASD")

# psdata n = 1500
psdata[1:3, ]

# SumStat Calculate summary statistics for propensity score weight
summary(msstat)
# importing user-supplied propensity scores "e.h"
fit <- nnet::multinom(formula=ps.formula, data=psdata, maxit=500, trace=FALSE)
e.h <- fit$fitted.values
varname <- c("cov1","cov2","cov3","cov4","cov5","cov6")
msstat0 <- SumStat(zname="trt", xname=varname, data=psdata, ps.estimate=e.h,
                   trtgrp="2", weight=c("ATE",'ATT',"ATO"))
summary(msstat0)


###########################
# PSweight Estimate causal effects by propensity score weighting

data("psdata")
# the propensity and outcome models
ps.formula <- trt~cov1+cov2+cov3+cov4+cov5+cov6
out.formula <- Y~cov1+cov2+cov3+cov4+cov5+cov6
# without augmentation
ato1<-PSweight(ps.formula = ps.formula,yname = 'Y',data = psdata,weight = 'ATO')
summary(ato1)
# augmented weighting estimator
ato2<-PSweight(ps.formula = ps.formula,yname = 'Y',data = psdata,
               augmentation = TRUE,out.formula = out.formula,family = 'gaussian',weight = 'ATO')
summary(ato2)

example("PSweight")


