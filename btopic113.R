## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE---------------------------------------
library(INLA)
rm(list=ls())
options(width=70, digits=2)
set.seed(2017)


## ------------------------------------------------------------------------
data(Seeds)
df = data.frame(y = Seeds$r, Ntrials = Seeds$n, Seeds[, 3:5])


## ------------------------------------------------------------------------
holdout = c(7, 12)
# - alternatively: sort(sample(1:nrow(df), 2))
df.holdout = df
df.holdout$y[holdout] = NA


## ------------------------------------------------------------------------
family1 = "binomial"
control.family1 = list(control.link=list(model="logit"))
# number of trials is df$Ntrials


## ------------------------------------------------------------------------
hyper1 = list(theta = list(prior="pc.prec", param=c(1,0.01)))
formula1 = y ~ x1 + x2 + f(plate, model="iid", hyper=hyper1)


## ------------------------------------------------------------------------
res1 = inla(formula=formula1, data=df.holdout, 
            family=family1, Ntrials=Ntrials, 
            control.family=control.family1,
            control.predictor=list(compute=T),
            control.compute=list(config=T))


## ------------------------------------------------------------------------
summary(res1)


## ------------------------------------------------------------------------
n.samples = 1000
samples = inla.posterior.sample(n.samples, result = res1)


## ------------------------------------------------------------------------
df[holdout, ]


## ------------------------------------------------------------------------
## See one example value
samples[[1]]$latent[holdout, , drop=F]

## Draw this part of every sample
samp = lapply(samples, function(x) x$latent[holdout])
samp = matrix(unlist(samp), ncol = 2, byrow = T)



## ------------------------------------------------------------------------
plot(density(samp[, 1]))
lines(density(samp[, 2]))


## ------------------------------------------------------------------------
res1$.args$control.family$control.link$model

samp.link = inla.link.invlogit(samp)


## ------------------------------------------------------------------------
plot(density(samp.link[, 1]))
## Add the naive estimate of binomial probability y/N:
abline(v = df$y[holdout[1]]/df$Ntrials[holdout[1]], col="blue")


## ------------------------------------------------------------------------
plot(density(samp.link[, 2]))
## Add the naive estimate of binomial probability y/N:
abline(v = df$y[holdout[2]]/df$Ntrials[holdout[2]], col="blue")


## ------------------------------------------------------------------------
## What range of observations do we want to compute
true.values = df$y[holdout]
if (any(is.na(true.values))) stop("Put the true values here!")

probs = rep(0, length(holdout))
for (i in 1:nrow(samp.link)) {
  probs = probs + dbinom(true.values, size=df$Ntrials[holdout], prob = samp.link[i, ])
}

## Numerical integration needs to divide by the number of samples
probs = probs / length(samples)
names(probs) = paste("holdout", holdout)

## Because of conditional independence, with the same sample, we can sum up the log probabilities
probs = c(probs, exp(sum(log(probs))))
names(probs)[length(probs)] = "Joint probability"

## The result for this hold-out combo
print(probs)

