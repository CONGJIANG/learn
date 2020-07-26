#############################################################################################
# Overlap weighting and sandwich variance with a logistic propensity score model
# Applicable to binary treatments, see the following Ref. for details

# Ref. Li F, Thomas LE, Li F. (2018) Addressing Extreme Propensity Scores via the 
# Overlap Weights. American Journal of Epidemiology

# This R program is intended to demonstrate the calculations of ATO and its variance
# as given in equations (2)-(4) of Li et al. (2018)

# 09/26/2018 The program has been extended to allow inference for causal risk ratio
#            and causal odds ratio with binary outcomes

# Our group is currently working on an R software that includes the functionality 
# of this demonstration program. 

# The following quantities are output from the program:
# (1) tau: estimator of the causal quantity among the overlap population
# (2) se:  empirical sandwich standard error
# (2) lcl: lower confidence limit
# (3) ucl: upper confidence limit
# (4) asd: absolute standardized difference (should be 0, see Web Appendix 2 of Li et al.)

# Frank Li; Sep 2018

# INPUT
# y: the n x 1 vector of outcomes
# z: the n x 1 vector of treatment indicators (e.g., trt = 1; control = 0)
# X: the n x p matrix of covariates for the PS model (excluding intercept term)
# estimand: "1" (ATO/risk difference), "2" (log risk ratio), "3" (log odds ratio)
#############################################################################################
OW <- function(y, z, X, estimand = 1){
  
  # module for checking balance in the weighted sample
  diff <- function(cov, z, ps){
    # cov: covariate
    # z: treatment status
    # ps: estimated propensity scores
    
    v1 <- cov[z == 1]
    v0 <- cov[z == 0]
    w1 <- 1-ps[z == 1]
    w0 <- ps[z == 0]
    n1 <- length(v1)
    n0 <- length(v0)
    delta <- abs(sum(v1*w1) / sum(w1) - sum(v0*w0) / sum(w0))
    # tstat <- delta / sqrt(var(v1)/n1 + var(v0)/n0)
    tstat <- delta / sqrt(var(v1) + var(v0))
    
    # return the absolute standardized difference
    return(tstat)
  }
  
  # summary statistics
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1             # total sample size
  W <- cbind(1, X)         # design matrix (including intercept)
  
  # estimate ps
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  
  # point estimate
  mu1.h <- sum(z*y*(1-e.h)) / sum(z*(1-e.h))
  mu0.h <- sum((1-z)*y*e.h) / sum((1-z)*e.h)
  
  # sandwich covariance estimator
  theta <- sum(e.h*(1-e.h)) / n
  H.b1 <- (z*(y-mu1.h)) * (e.h*(1-e.h)) * W
  H.b0 <- (-(1-z)*(y-mu0.h)) * (e.h*(1-e.h)) * W
  H.b1 <- colSums(H.b1) / n
  H.b0 <- colSums(H.b0) / n
  
  # Note: the following information matrix can also be extracted 
  # directly from the 'fit' object. Here we compute it from first 
  # principles
  E.bb <- crossprod(sqrt(e.h*(1-e.h)) * W) / n   
  psi1 <- z*(1-e.h)*(y-mu1.h) - 
    (z-e.h) * c(t(H.b1) %*% solve(E.bb) %*% t(W))
  psi0 <- (1-z)*e.h*(y-mu0.h) - 
    (z-e.h) * c(t(H.b0) %*% solve(E.bb) %*% t(W))
  Psi <- rbind(psi1, psi0)
  V <- tcrossprod(Psi) / (n*theta)^2
  
  # Compute causal quantity
  # we use the delta method for log RR and log OR to improve 
  # the normality approximation (this is what most software does, 
  # as RR and OR tends to be skewed)
  # one could simply exponentiate the confidence limits to get
  # confidence limits for RR and OR
  if(estimand == 1){                     # ATO/risk diff
    tau <- mu1.h - mu0.h
    a <- c(1,-1)
  } else if(estimand == 2){              # log risk ratio
    tau <- log(mu1.h) - log(mu0.h)
    a <- c(1/mu1.h, -1/mu0.h) 
  } else{                                # log odds ratio
    tau <- log(mu1.h/(1-mu1.h)) - log(mu0.h/(1-mu0.h))
    a <- c(1/mu1.h + 1/(1-mu1.h), -1/mu0.h - 1/(1-mu0.h)) 
  }
  vtau10 <- t(a) %*% V %*% a
  se <- sqrt(vtau10)
  # Confidence limits
  lcl <- tau - qnorm(0.975)*se
  ucl <- tau + qnorm(0.975)*se
  
  # check balance (expected to be zero)
  asd <- rep(NA, ncol(X))
  for(k in 1:ncol(X)){
    asd[k] <- diff(cov = X[,k], z = z, ps = e.h)
  }
  
  # output the quantities of interest
  return(list(tau=tau, se=se, lcl=lcl, ucl = ucl, asd = asd))
}

###########################################################
# Illustrative calculations (1) - ATO
###########################################################
# simulate a data set with six potential confounders
require(mvtnorm)

n <- 1000
set.seed(123)
X <- rmvnorm(n, rep(0,6), diag(6))
W <- cbind(1, X)

# true propensity scores
e <- plogis(c(W%*%c(0.3, rep(0.25,6))))
z <- rbinom(n, 1, e)

# generate outcomes (assume additive for simplicity)
delta <- 1
mu.y <- c(W%*%c(1,rep(0.5,3),rep(-0.5,3))) + z*delta
y <- rnorm(n, mu.y, sd = 2)

# use the function
res <- OW(y=y, z=z, X=X, estimand=1)
res$tau
res$se

# compare with canonical treatment model
require(survey)

ehat <- as.numeric(glm(z ~ X,family=binomial(link="logit"))$fitted.values)
w <- rep(NA, n)
w[z == 1] <- 1-ehat[z == 1]
w[z == 0] <- ehat[z == 0]
w[z == 1] <- w[z == 1]/sum(w[z == 1])
w[z == 0] <- w[z == 0]/sum(w[z == 0])
data <- data.frame(y=y, z=z)
design <- svydesign(ids=~1, weights=~w, data=data)

tab <- summary(svyglm(y ~ z,design=design,family=gaussian(link="identity")))$coef
tab[2,1]
tab[2,2]  # conservative variance

###########################################################
# Illustrative calculations (2) - ratio
###########################################################
n <- 1000
set.seed(123)
X <- rmvnorm(n, rep(0,6), diag(6))
W <- cbind(1, X)

# true propensity scores
e <- plogis(c(W%*%c(0.3, rep(0.25,6))))
z <- rbinom(n, 1, e)

# generate outcomes (assume additive for simplicity)
delta <- 1
p1 <- plogis(c(W%*%c(1,rep(1,3),rep(-0.5,3))))
p0 <- plogis(c(W%*%c(1,rep(0,3),rep(0.1,3))))
y1 <- rbinom(n, 1, p1)
y0 <- rbinom(n, 1, p0)
y <- y1
y[z == 0] <- y0[z == 0]

# true value
TAU1 <- mean(e*(1-e)*p1)/mean(e*(1-e))
TAU0 <- mean(e*(1-e)*p0)/mean(e*(1-e))
log(TAU1/TAU0)
log(TAU1/(1-TAU1)) - log(TAU0/(1-TAU0))

# use the function
res_lRR = OW(y=y, z=z, X=X, estimand=2)
res_lOR = OW(y=y, z=z, X=X, estimand=3)

res_lRR$tau
res_lRR$se

res_lOR$tau
res_lOR$se

# compare with canonical treatment model
ehat <- as.numeric(glm(z ~ X,family=binomial(link="logit"))$fitted.values)
w <- rep(NA, n)
w[z == 1] <- 1-ehat[z == 1]
w[z == 0] <- ehat[z == 0]
w[z == 1] <- w[z == 1]/sum(w[z == 1])
w[z == 0] <- w[z == 0]/sum(w[z == 0])
data <- data.frame(y=y, z=z)
design <- svydesign(ids=~1, weights=~w, data=data)

tab_lRR <- summary(svyglm(y ~ z,design=design,family=binomial(link="log")))$coef
tab_lRR[2,1]
tab_lRR[2,2]  # conservative variance

tab_lOR <- summary(svyglm(y ~ z,design=design,family=binomial(link="logit")))$coef
tab_lOR[2,1]
tab_lOR[2,2]  # conservative variance