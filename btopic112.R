## ----setup, include=FALSE-----------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE--------------------------------
library(INLA)
rm(list=ls())
options(width=70, digits=2)
set.seed(2017)


## -----------------------------------------------------------------
data(Seeds)
df = data.frame(y = Seeds$r, Ntrials = Seeds$n, Seeds[, 3:5])


## -----------------------------------------------------------------
family1 = "binomial"
control.family1 = list(control.link=list(model="logit"))
# number of trials is df$Ntrials


## -----------------------------------------------------------------
hyper1 = list(theta = list(prior="pc.prec", param=c(1,0.01)))
formula1 = y ~ x1 + x2 + f(plate, model="iid", hyper=hyper1)


## -----------------------------------------------------------------
res1 = inla(formula=formula1, data=df, 
            family=family1, Ntrials=Ntrials, 
            control.family=control.family1,
            control.predictor=list(compute=T),
            control.compute=list(config=T))


## -----------------------------------------------------------------
summary(res1)


## -----------------------------------------------------------------
n.samples = 1000
samples = inla.posterior.sample(n.samples, result = res1)


## -----------------------------------------------------------------
str(samples[[1]])


## -----------------------------------------------------------------
t(samples[[1]]$latent)
# - transposed for shorter printing
# - this is a column vector


## -----------------------------------------------------------------
res1$misc$configs$contents


## -----------------------------------------------------------------
contents = res1$misc$configs$contents
effect = "plate"
id.effect = which(contents$tag==effect)
# - the numerical id of the effect
ind.effect = contents$start[id.effect]-1 + (1:contents$length[id.effect])
# - all the indices for the effect
# - these are the indexes in the sample[[1]]$latent !


## -----------------------------------------------------------------
## See an example for the first sample
samples[[1]]$latent[ind.effect, , drop=F]

## Draw this part of every sample
samples.effect = lapply(samples, function(x) x$latent[ind.effect])


## -----------------------------------------------------------------
s.eff = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
# - s.eff means samples.effect, just as a matrix
colnames(s.eff) = rownames(samples[[1]]$latent)[ind.effect]
# - retrieve names from original object
summary(s.eff)


## -----------------------------------------------------------------
cbind(colMeans(s.eff), res1$summary.random$plate$mean)


## -----------------------------------------------------------------
df[7, ]


## -----------------------------------------------------------------
## see one example value
samples[[1]]$latent[7, , drop=F]

## Draw this part of every sample
samples.pred7 = lapply(samples, function(x) x$latent[7])
samples.pred7 = unlist(samples.pred7)



## -----------------------------------------------------------------
plot(density(samples.pred7))
abline(v=df$y[7])


## -----------------------------------------------------------------
res1$.args$control.family$control.link$model


## -----------------------------------------------------------------
samples.link7 = inla.link.invlogit(samples.pred7)

plot(density(samples.link7))

## Add the naive estimate of binomial probability y/N:
abline(v = df$y[7]/df$Ntrials[7], col="blue")



## -----------------------------------------------------------------
## What range of observations do we want to compute
discrete.range = 0:100
## Initialize the probabilities 
probs = rep(0, length(discrete.range))

for (i in 1:length(samples.link7)) {
  probs = probs + dbinom(discrete.range, size=df$Ntrials[7], prob = samples.link7[i])
}
probs = probs / length(samples)
names(probs) = discrete.range


## -----------------------------------------------------------------
plot(probs, xlim=c(37, 67))
abline(v=df$y[7], col="blue")


## -----------------------------------------------------------------
## Remember the covariate values
df[7, ]

## Remember the formula
res1$.args$formula


## -----------------------------------------------------------------
nr = 57
s = samples[[nr]]$latent

## beta1 * 0
f.x1.0 = s[44, , drop=F] * 0
## beta1 * 1
f.x2.1 = s[45, , drop=F] * 1
## f(plate = 7)
f.plate.7 = s[28, , drop=F]
## The intercept
int = s[43, , drop=F]

sum = drop(f.x1.0 + f.x2.1 + f.plate.7 + int)
list(model.component.sum.7 = sum, linear.predictor.7 = samples.pred7[nr])



## -----------------------------------------------------------------
our.experiment = list(x1 = 0, x2 = 0, plate = 7)

nr = 57
s = samples[[nr]]$latent

## beta1 * 0
f.x1.0 = s[44, , drop=F] * our.experiment$x1
## beta1 * 1
f.x2.1 = s[45, , drop=F] * our.experiment$x2
## f(plate = 7)
# - the same plate as before
f.plate.7 = s[28, , drop=F]
## The intercept
int = s[43, , drop=F]

predictor.our.experiment = drop(f.x1.0 + f.x2.1 + f.plate.7 + int)


## -----------------------------------------------------------------
## As before:
our.experiment = list(x1 = 0, x2 = 0, plate = 50)
nr = 57

## Want a hyper-parameter
plate.precision = samples[[nr]]$hyperpar[1]
plate.sigma = plate.precision^-0.5
names(plate.sigma) = "sigma"

print(plate.sigma)



## -----------------------------------------------------------------
sample.plate50 = rnorm(1, sd=plate.sigma)


## -----------------------------------------------------------------
## As before
s = samples[[nr]]$latent
f.x1.0 = s[44, , drop=F] * our.experiment$x1
f.x2.1 = s[45, , drop=F] * our.experiment$x2

## f(plate = 50)
f.plate.50 = sample.plate50
## The intercept
int = s[43, , drop=F]

predictor.experiment.newplate = drop(f.x1.0 + f.x2.1 + f.plate.50 + int)

