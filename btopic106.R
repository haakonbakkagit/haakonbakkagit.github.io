## ----setup, include=FALSE-----------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE--------------------------------
library(INLA); rm(list=ls())
options(width=70, digits=2)


## -----------------------------------------------------------------
data(Seeds); dat = Seeds
df = data.frame(y = dat$r, Ntrials = dat$n, dat[, 3:5])


## -----------------------------------------------------------------
family1 = "binomial"
control.family1 = list(control.link=list(model="logit"))


## -----------------------------------------------------------------
hyper1 = list(theta = list(prior="pc.prec", param=c(1,0.01)))
formula1 = y ~ f(plate, model="iid", hyper=hyper1)


## -----------------------------------------------------------------
res1 = inla(formula=formula1, data=df, 
            family=family1, Ntrials=Ntrials, 
            control.family=control.family1,
            control.compute=list(config=T))
# - use control.compute so that we get to sample from the posterior


## -----------------------------------------------------------------
summary(res1)


## -----------------------------------------------------------------
(avg.all = mean(res1$summary.random$plate$mean))
(avg.x1.0 = mean(res1$summary.random$plate$mean[df$x1==0]))
(avg.x1.1 = mean(res1$summary.random$plate$mean[df$x1==1]))
(avg.diff = avg.x1.0 - avg.x1.1)


## -----------------------------------------------------------------
n.samples = 10000
samples = inla.posterior.sample(n.samples, result = res1)
(mean(samples[[1]]$latent[(1:nrow(df))][df$x1==0]) - mean(samples[[1]]$latent[(1:nrow(df))][df$x1==1]))
# - this gives the average difference between the linear predictor values for plates with the two different seed types
samples.avg.diff = unlist(lapply(samples, FUN = function(x) mean(x$latent[(1:nrow(df))][df$x1==0]) - mean(x$latent[(1:nrow(df))][df$x1==1])))
# - this uses the line above in lapply, to get all the samples
summary(samples.avg.diff)
# - Note that the mean equals avg.diff


## -----------------------------------------------------------------
plot(density(samples.avg.diff))

