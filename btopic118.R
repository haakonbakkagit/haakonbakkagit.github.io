## ----setup, include=FALSE--------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- message=FALSE, warning=FALSE-----------------------------
library(INLA)
set.seed(2017)


## --------------------------------------------------------------
if (T) {
  n.group = 10
} else {
  n.group = 40
}
group = rep(1:n.group, times=(1:n.group))
# - 1 individual in group 1, 2 in group 2 etc
sig.group = 2
sig.ind = 0.5
group.effect = drop(scale(1:n.group))*sig.group
group.effect.i = group.effect[group]
N = length(group)
individual.effect = rnorm(N)*sig.ind

eta = 0.5+individual.effect + group.effect.i
y = rpois(n=N, lambda = exp(eta))


## --------------------------------------------------------------
df = data.frame(y = y, i = 1:N, group = group)


## --------------------------------------------------------------
hyper.iid = list(theta = list(prior="pc.prec", param=c(3, 0.5)))
# - this is a very weak prior

formula = y ~ f(i, model="iid", hyper=hyper.iid)+ f(group, model="iid", hyper=hyper.iid)

res = inla(formula, data=df, family="Poisson", E = rep(1, N),
           control.predictor=list(compute=TRUE, link=1))


## --------------------------------------------------------------
summary(res)


## --------------------------------------------------------------
plot(df$y, col="blue", pch=1)
points(res$summary.fitted.values$mean, pch=8)
points(exp(eta), pch=20, col="blue")
for (j in 1:n.group) {
  abline(v = which(df$group==j)[1]-0.5, col="grey")
}


## --------------------------------------------------------------
plot(res$summary.fitted.values$mean, df$y)

