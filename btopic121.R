## ----setup, include=FALSE--------------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- message=FALSE, warning=FALSE-----------------------------------------------
library(INLA)


## --------------------------------------------------------------------------------
set.seed(2017)
n = 50
idx = 1:n
fun = 100*((idx-n/2)/n)^3
y = fun + rnorm(n, mean=0, sd=1)
# - add some noise
plot(idx, y)
lines(fun, col="darkgreen")
df = data.frame(y=y, idx=idx)


## --------------------------------------------------------------------------------
hyper1 = list(prec=list(initial=0, fixed=T))
formula1 = y ~ -1 + f(idx, model="rw2", hyper=hyper1)
res1 = inla(formula1, family="gaussian", data=df,
            control.family = list(hyper = hyper1))


## --------------------------------------------------------------------------------
local.plot.result = function(res) {
  plot(idx, y)
  lines(res$summary.random$idx$mean, col="blue")
  lines(res$summary.random$idx$`0.025quant`, col="grey")
  lines(res$summary.random$idx$`0.9`, col="grey")
  lines(fun, col="darkgreen")
}
local.plot.result(res1)


## --------------------------------------------------------------------------------
formula2 = y ~ -1 + f(idx, model="rw2")
res2 = inla(formula2, family="gaussian", data=df,
            control.family = list(hyper = hyper1),
            control.compute = list(config=TRUE))


## --------------------------------------------------------------------------------
plot(res2$marginals.hyperpar$`Precision for idx`, type="l", xlim=c(0, 120))


## --------------------------------------------------------------------------------
plot(res2$internal.marginals.hyperpar$`Log precision for idx`, type="l")


## --------------------------------------------------------------------------------
str(res2$misc$configs$config[[1]], 1)


## --------------------------------------------------------------------------------
data.frame( theta = unlist(lapply(res2$misc$configs$config, function(x) x$theta)),
       log.post = unlist(lapply(res2$misc$configs$config, function(x) x$log.posterior)))


## --------------------------------------------------------------------------------
local.plot.result(res2)


## --------------------------------------------------------------------------------
res3 = inla(formula2, family="gaussian", data=df)


## --------------------------------------------------------------------------------
local.plot.result(res2)

