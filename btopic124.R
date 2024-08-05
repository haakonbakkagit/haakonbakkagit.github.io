## ---- include=FALSE-----------------------------------------------
rm(list=ls())
options(width=70, digits=2)
## TODOS
# (ignore these comments)
#For an overview of all the related topics see [TODO].
#For an overview of the words and concepts see [TODO].


## ---- warning=FALSE, message=FALSE--------------------------------
library(INLA)
library(fields)
library(ggplot2)
set.seed(201803)


## -----------------------------------------------------------------
N = 500
x1 = runif(N, min=-2, max=2)
x2 = runif(N, min=-2, max=2)
x3 = runif(N, min=-2, max=2)
beta0 = 3
beta1 = 0.7
beta2 = -0.3
v1 = 0.3
sig.eps = 0.4


## -----------------------------------------------------------------
eta = beta0 + beta1*x1 + beta2*x2 + v1*x3
y = eta + sig.eps*rnorm(N)


## -----------------------------------------------------------------
df = data.frame(y = y, x1=x1, x2=x2)


## -----------------------------------------------------------------
formula = y ~ x1 + x2
res1 = inla(formula, family = "gaussian", data = df,
            control.predictor = list(compute=T))


## -----------------------------------------------------------------
res1$summary.fixed


## -----------------------------------------------------------------
(sd.model.noise = sd(y-res1$summary.linear.predictor$mean))
(sd.true.noise = sd(y-eta))
# - approximately equal to sig.eps
(sd.v1x3 = sd(v1*x3))
(sqrt(sd.true.noise^2 + sd.v1x3^2))


## -----------------------------------------------------------------
formula = y ~ x1 + x2 + x3
res2 = inla(formula, family = "gaussian", data = df,
            control.predictor = list(compute=T))


## -----------------------------------------------------------------
res2$summary.fixed

