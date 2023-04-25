## ----setup, include=FALSE--------------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## --------------------------------------------------------------------------------
## Empirical quantile transform
x = c(0:10, 20:21, 100:105, 1000)
qtx = ecdf(x)
## The transformed variable:
qtx(x)


## ---- message=FALSE, warning=FALSE-----------------------------------------------
library(INLA)
inla.setOption("num.threads", 1)
set.seed(2021)


## --------------------------------------------------------------------------------
## Number of observations
n.obs = 5E3
## The Gaussian noise
sig.epsilon = 1.00
## Covariates
n.cov.important = 7
n.cov.noeffect = 20
## Fraction of new variable that is taken from some previous variable
## 0 gives uncorrelated, infinity gives perfectly correlated
corr.fraction = 0.7
## How many factors
## Please do not change this number (that requires further code changes)
n.factors = 6
## How many factor levels, at most
n.fact.levels = 7
## Typical effect size of important factor level
effect.size.fact = 2.10


## --------------------------------------------------------------------------------
## Matrix of covariates 
## Build this column by column
X = cbind(rnorm(n.obs))
## What previous variable are you correlated with?
corr.with = rep(NA, n.cov.important)
for (i in 2:(n.cov.important+n.cov.noeffect)) {
  ## Choose a previous covariate to be correlated with
  corr.with[i] = sample(1:(i-1), 1)
  ## Create a new covariate
  new = (corr.fraction*X[, corr.with[i]] + rnorm(n.obs))/corr.fraction
  ## Add a column to the matrix of covariates
  X = cbind(X, new)
}
colnames(X) = paste0("x", 1:ncol(X))


## --------------------------------------------------------------------------------
summary(X[, 1:5])


## --------------------------------------------------------------------------------
cor(X[, 1:5])


## --------------------------------------------------------------------------------
## Simulate a factor with n.fact.levels factor levels
new = sample(1:n.fact.levels, size=n.obs, replace = T)
## The factor (as strings)
new = paste0("F", 1, "level", new)
## Create the first column of the matrix of factors
FF = cbind(new)
## The effect of the factor
FFeff = cbind(rep(0, n.obs))
## Each factor has an effect level and an effect size
## Effect level gives which of the factor levels has a nonzero effect
effect.level = rep(NA, n.factors)
## Effect size gives how large the effect for the nonzero factor is
## The first factor has an effect of 0
effect.size = rep(NA, n.factors)
for (i in 2:(n.factors)) {
  ## Choose how many factor levels for this factor
  tmp.levels = sample(2:n.fact.levels, 1)
  ## Create a new factor
  new = sample(1:tmp.levels, size=n.obs, replace = T)
  ## Create its effect
  effect.level[i] = sample(1:tmp.levels, 1)
  effect.size[i] = rnorm(1)*effect.size.fact
  new.effect = (new==effect.level[i])*effect.size[i]
  ## Create the column for the new factor and add it to the matrix
  new = paste0("F", i, "level", new)
  FF = cbind(FF, new)
  ## Add the new effect column to the matrix
  FFeff = cbind(FFeff, new.effect)
}
colnames(FF) = paste0("fa", 1:ncol(FF))
colnames(FFeff) = paste0("fa", 1:ncol(FF), "eff")


## --------------------------------------------------------------------------------
FF[1:8, 2:6]

## --------------------------------------------------------------------------------
FFeff[1:8, 2:6]


## --------------------------------------------------------------------------------
## 2 alternatives, no effect at all or negligble effect
betas = c(abs(rnorm(n.cov.important))+0.5, rnorm(n.cov.noeffect)*0.000) * 0.1
# betas = c(abs(rnorm(n.cov.important))+0.5, rnorm(n.cov.noeffect)*0.001) * 0.1


## --------------------------------------------------------------------------------
gaussian.error = rnorm(n.obs, 0, sd=sig.epsilon)


## --------------------------------------------------------------------------------
## Compute the response from the simulated values
eta1 = X %*% betas
eta2 = rowSums(FFeff)
y = eta1 + eta2 + gaussian.error


## --------------------------------------------------------------------------------
sd(eta1)
sd(eta2)
sd(gaussian.error)


## --------------------------------------------------------------------------------
form1 = y ~ 1 + X
fit1 = inla(form1, family="gaussian", data = list(y=y, X=X))


## --------------------------------------------------------------------------------
## Look at the first 10 fixed effects
fit1$summary.fixed[1:8, 1:2]

## --------------------------------------------------------------------------------
## How far is the estimate from the truth?
betas - fit1$summary.fixed$mean[-1]


## --------------------------------------------------------------------------------
## Total RMSE
sqrt(sum((betas - fit1$summary.fixed$mean[-1])^2))
## Is the estimate within two std errors?
table((abs(betas - fit1$summary.fixed$mean[-1]))<fit1$summary.fixed$sd[-1]*2)
## Quantiled model sd
quantile(fit1$summary.fixed$sd[-1])


## --------------------------------------------------------------------------------
fit2 = inla(form1, family="gaussian", data = list(y=y, X=X),
            control.fixed = list(prec=1))


## --------------------------------------------------------------------------------
## Total RMSE
sqrt(sum((betas - fit2$summary.fixed$mean[-1])^2))
## Is the estimate within two std errors?
table((abs(betas - fit2$summary.fixed$mean[-1]))<fit2$summary.fixed$sd[-1]*2)
## Quantiled model sd
quantile(fit2$summary.fixed$sd[-1])


## --------------------------------------------------------------------------------
fit3 = inla(form1, family="gaussian", data = list(y=y, X=X),
            control.fixed = list(prec=100))


## --------------------------------------------------------------------------------
## Total RMSE
sqrt(sum((betas - fit3$summary.fixed$mean[-1])^2))
## Is the estimate within two std errors?
table((abs(betas - fit3$summary.fixed$mean[-1]))<fit3$summary.fixed$sd[-1]*2)
## Quantiled model sd
quantile(fit3$summary.fixed$sd[-1])


## --------------------------------------------------------------------------------
## Fixed hyperparameter in the iid effect is the same as putting
## control.fixed = list(prec="what you fix it to")
## The theta1 is the log precision
hyper.fix = list(theta1 = list(initial=log(1), fixed=T))
## Dynamic hyperparameter, with a prior for log precision that has median 0.1
hyper.iid = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))


## --------------------------------------------------------------------------------
data1 = as.list(as.data.frame(FF))
data1$y = drop(y)
str(data1)


## --------------------------------------------------------------------------------
form2 = y ~ 1 + f(fa1, model="iid", hyper=hyper.fix) + f(fa2, model="iid", hyper=hyper.fix) + f(fa3, model="iid", hyper=hyper.fix) + f(fa4, model="iid", hyper=hyper.fix) + f(fa5, model="iid", hyper=hyper.fix) + f(fa6, model="iid", hyper=hyper.fix)
fit4 = inla(form2, family="gaussian", data = data1)


## --------------------------------------------------------------------------------
fit4$summary.random$fa1[, 1:3]


## --------------------------------------------------------------------------------
print(paste("Effect on level", effect.level[2], "is", effect.size[2], 
            "compared to baseline 0"))
fit4$summary.random$fa2[, 1:3]


## --------------------------------------------------------------------------------
form3 = y ~ 1 + f(fa1, model="iid", hyper=hyper.iid) + f(fa2, model="iid", hyper=hyper.iid) + f(fa3, model="iid", hyper=hyper.iid) + f(fa4, model="iid", hyper=hyper.iid) + f(fa5, model="iid", hyper=hyper.iid) + f(fa6, model="iid", hyper=hyper.iid)
fit5 = inla(form3, family="gaussian", data = data1,
            control.inla = list(int.strategy="eb"))


## --------------------------------------------------------------------------------
## We see that the precision for fa1 is very big, so there is little effect of fa1
fit5$summary.hyperpar[, 1:2]


## --------------------------------------------------------------------------------
## The following uncertainties are much better than in model 4
fit5$summary.random$fa1[, 1:3]


## --------------------------------------------------------------------------------
## The following uncertainties are better than in model 4
fit5$summary.random$fa2[, 1:3]


## --------------------------------------------------------------------------------
form4 = y ~ 1 + X + f(fa1, model="iid", hyper=hyper.iid) + f(fa2, model="iid", hyper=hyper.iid) + f(fa3, model="iid", hyper=hyper.iid) + f(fa4, model="iid", hyper=hyper.iid) + f(fa5, model="iid", hyper=hyper.iid) + f(fa6, model="iid", hyper=hyper.iid)
fit6 = inla(form4, family="gaussian", data = c(data1, list(X=X)),
            control.inla = list(int.strategy="eb"),
            control.fixed = list(prec=100))


## --------------------------------------------------------------------------------
## Total RMSE
sqrt(sum((betas - fit6$summary.fixed$mean[-1])^2))
## Is the estimate within two std errors?
table((abs(betas - fit6$summary.fixed$mean[-1]))<fit6$summary.fixed$sd[-1]*2)
## Quantiled model sd
quantile(fit6$summary.fixed$sd[-1])


## --------------------------------------------------------------------------------
print(paste("Effect on level", effect.level[2], "is", effect.size[2], 
            "compared to baseline 0"))
fit6$summary.random$fa2[, 1:3]

