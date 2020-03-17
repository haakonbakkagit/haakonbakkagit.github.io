## ---- include=FALSE-------------------------------------------------
rm(list=ls())
options(width=70, digits=2)
## TODOS
# (ignore these comments)
#  


## ---- warning=FALSE, message=FALSE----------------------------------
library(INLA)
library(fields)
library(ggplot2)
set.seed(201803)


## -------------------------------------------------------------------
N = 200
x1 = runif(N, min=-2, max=2)
x2 = runif(N, min=-2, max=2)
x3 = runif(N, min=-2, max=2)
beta0 = 3
beta1 = 0.7
beta2 = -0.3
v1 = 0.3


## -------------------------------------------------------------------
eta = beta0 + beta1*x1 + beta2*x2 + v1*x3
y = rpois(n=N, lambda = exp(eta))


## -------------------------------------------------------------------
df = data.frame(y = y, x1=x1, x2=x2)


## -------------------------------------------------------------------
formula1 = y ~ x1 + x2
res1 = inla(formula1, family = "poisson", data = df,
            control.predictor = list(compute=T))


## -------------------------------------------------------------------
res1$summary.fixed


## -------------------------------------------------------------------
formula2 = y ~ x1 + x2 + f(id.iid, model="iid")
res2 = inla(formula2, family = "poisson", 
            data = data.frame(df, id.iid = 1:nrow(df)),
            control.predictor = list(compute=T))


## -------------------------------------------------------------------
res2$summary.fixed


## -------------------------------------------------------------------
formula3 = y ~ x1 + x2 + x3
res3 = inla(formula3, family = "poisson", 
            data = data.frame(df, id.iid = 1:nrow(df)),
            control.predictor = list(compute=T))


## -------------------------------------------------------------------
res3$summary.fixed


## -------------------------------------------------------------------
one.sim.coverage = function(x) {
  ## sim data
  df$y = rpois(n=N, lambda = exp(eta))
  
  ## inference
  res1 = inla(formula1, family = "poisson", data = df,
              quantiles = c(0.1, 0.9),
              num.threads = 1)
  res2 = inla(formula2, family = "poisson", 
              data = data.frame(df, id.iid = 1:nrow(df)),
              quantiles = c(0.1, 0.9), 
              num.threads = 1)
  
  ## coverage T/F of 80% interval
  a1 = (beta1 > res1$summary.fixed$`0.1quant`[2] & beta1 < res1$summary.fixed$`0.9quant`[2])
  a1.2 = (beta2 > res1$summary.fixed$`0.1quant`[3] & beta2 < res1$summary.fixed$`0.9quant`[3])
  a2 = (beta1 > res2$summary.fixed$`0.1quant`[2] & beta1 < res2$summary.fixed$`0.9quant`[2])
  a2.2 = (beta2 > res2$summary.fixed$`0.1quant`[3] & beta2 < res2$summary.fixed$`0.9quant`[3])
  return(c(a1, a1.2, a2, a2.2))
}


## -------------------------------------------------------------------
n.sim = 100
if (require(parallel)) {
  # - run in parallel if you have this package installed
  res.sim.raw = mclapply(as.list(1:n.sim), mc.cores=3, FUN = one.sim.coverage)
} else {
  res.sim.raw = lapply(as.list(1:n.sim), FUN = one.sim.coverage)
}


## -------------------------------------------------------------------
res.sim = as.matrix(data.frame(res.sim.raw))
colnames(res.sim) = NULL
## Model 1: beta1
table(res.sim[1, ])
## Model 1: beta2
table(res.sim[2, ])
## Model 2: beta1
table(res.sim[3, ])
## Model 2: beta2
table(res.sim[4, ])

