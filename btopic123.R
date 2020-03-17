## ----setup, include=FALSE-------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- include=FALSE-------------------------------------------------
rm(list=ls())
options(width=70, digits=2)
## TODOS
# (ignore these comments)
#For an overview of all the related topics see [TODO].
#For an overview of the words and concepts see [TODO].


## ---- warning=FALSE, message=FALSE----------------------------------
library(INLA)


## -------------------------------------------------------------------
data(Seeds)


## -------------------------------------------------------------------
# Run: ?Seeds
head(Seeds)
# - r is the number of seed germinated (successes)
# - n is the number of seeds attempted (trials)
# - x1 is the type of seed
# - x2 is the type of root extract
# - plate is the numbering of the plates/experiments


## -------------------------------------------------------------------
df = data.frame(y = Seeds$r, Ntrials = Seeds$n, Seeds[, 3:5])


## -------------------------------------------------------------------
family1 = "binomial"
control.family1 = list(control.link=list(model="logit"))
# number of trials is df$Ntrials


## -------------------------------------------------------------------
hyper1 = list(theta = list(prior="pc.prec", param=c(1,0.01)))
formula1 = y ~ x1 + x2 + f(plate, model="iid", hyper=hyper1)


## -------------------------------------------------------------------
res1 = inla(formula=formula1, data=df, 
            family=family1, Ntrials=Ntrials, 
            control.family=control.family1)


## -------------------------------------------------------------------
m.sigma = inla.tmarginal(fun = function(x) exp(-1/2*x), marginal = 
                           res1$internal.marginals.hyperpar$`Log precision for plate`)
# - m.sigma is the marginal for the standard deviation parameter in the iid random effect
plot(m.sigma, type="l", xlab = expression(sigma[iid]), ylab = "Probability density")


## -------------------------------------------------------------------
m.sigma = inla.tmarginal(fun = function(x) exp(x), marginal = 
                           res1$internal.marginals.hyperpar$`Log precision for plate`)
plot(m.sigma, type="l", xlab = "precision", ylab = "Probability density", xlim=c(0, 150))


## -------------------------------------------------------------------
m.sigma = inla.tmarginal(fun = function(x) exp(-x), marginal = 
                           res1$internal.marginals.hyperpar$`Log precision for plate`)
plot(m.sigma, type="l", xlab = "variance", ylab = "Probability density")


## -------------------------------------------------------------------
fun = function(x) -1*atan((3-10*x)*1.2)
xval = -10000:10000/10000 
plot(xval, fun(xval))


## -------------------------------------------------------------------
m.sigma = inla.tmarginal(fun = function(x) exp(-1/2*x), marginal = 
                           res1$internal.marginals.hyperpar$`Log precision for plate`)
m.sigma = inla.tmarginal(fun = fun, marginal = m.sigma)
plot(m.sigma, type="l", xlab = "special 1", ylab = "Probability density")


## -------------------------------------------------------------------
## The derivative
# - always positive
fun = function(x) sin(31*(x+0.55)) + 1.1
xval = 0:10000/10000 
plot(xval, fun(xval))

## The transformation function
# - Integral of the previous
# - Monotone
fun = function(x) -1/31*cos(31*(x+0.55)) + 1.1*x
xval = 0:10000/10000 
plot(xval, fun(xval))


## -------------------------------------------------------------------
m.sigma = inla.tmarginal(fun = function(x) exp(-1/2*x), marginal = 
                           res1$internal.marginals.hyperpar$`Log precision for plate`)
m.sigma = inla.tmarginal(fun = fun, marginal = m.sigma)
plot(m.sigma, type="l", xlab = "Special 2", ylab = "Probability density")

