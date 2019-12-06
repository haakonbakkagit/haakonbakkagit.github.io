## ----setup, include=FALSE------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- message=FALSE, warning=FALSE---------------------------------------
library(INLA)
inla.setOption("num.threads", 1)
# - to only run 1 processor thread (no paralellisation)
set.seed(2017)


## ------------------------------------------------------------------------
N = 3600
# - the number of observations
# - should be a multiple of 12
sig.epsilon = 0.5
# - the Gaussian noise
sig.u = 1
# - the structured part of the time series
rho = 0.90
# - the true autocorrelation
sig.seasonal = 2
# - the size of the seasonal effect


## ------------------------------------------------------------------------
u = arima.sim(list(order = c(1,0,0), ar = rho), n = N,sd=1)
# - this sd is not the marginal standard deviation
u = u/sd(u)*sig.u
# - this has the correct standard deviation


## ------------------------------------------------------------------------
seas.coeff = (0:11)*(1:12-12)
seas = rep(seas.coeff, N/12)
seas = drop(scale(seas))*sig.seasonal


## ------------------------------------------------------------------------
gaussian.error = rnorm(N, 0, sd=sig.epsilon)


## ------------------------------------------------------------------------
y = u + seas + gaussian.error


## ------------------------------------------------------------------------
df = data.frame(y = y, t = 1:N, year = rep(1:12, N/12))


## ------------------------------------------------------------------------
plot(df$t, df$y, main="Data", col="blue")


## ------------------------------------------------------------------------
hyper.ar1 = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)),
                 theta2 = list(prior="pc.cor1", param=c(0.9, 0.5)))
hyper.rw2 = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))
hyper.family = list(theta = list(prior="pc.prec", param=c(3, 0.5)))


## ------------------------------------------------------------------------
formula <- y~ f(t,model='ar1', hyper=hyper.ar1) + f(year, model="rw2", hyper=hyper.rw2, cyclic=T, constr=T)

res <- inla(formula=formula, data=df, family="gaussian",
             control.predictor=list(compute=TRUE),
             control.family = list(hyper = hyper.family))


## ------------------------------------------------------------------------
summary(res)


## ------------------------------------------------------------------------
plot(df$y, res$summary.fitted.values$mean, main="Fitting result")


## ------------------------------------------------------------------------
marginal = res$internal.marginals.hyperpar$`Log precision for the Gaussian observations`
transform = function(x) exp(-0.5*x)
sig.eps.posterior = inla.tmarginal(transform, marginal)

plot(sig.eps.posterior, type="l", xlab = expression(sigma), ylab = "Probability density",
     main = "Size of noise component")
    
xvals = seq(0.45, 0.62, length.out=1000)
lambda = -log(hyper.family$theta$param[2])/hyper.family$theta$param[1]
lines(xvals, 1E1*exp(-lambda*xvals), lty='dashed')
abline(v=sig.epsilon, col="blue")


## ------------------------------------------------------------------------
marginal = res$internal.marginals.hyperpar$`Log precision for t`
transform = function(x) exp(-0.5*x)
sig.posterior = inla.tmarginal(transform, marginal)

plot(sig.posterior, type="l", xlab = expression(sigma), ylab = "Probability density",
     main = "Size of AR1 component")
    
xvals = seq(0.5, 1.5, length.out=1000)
lambda = -log(hyper.ar1$theta1$param[2])/hyper.ar1$theta1$param[1]
lines(xvals, 1E3*exp(-lambda*xvals), lty='dashed')
abline(v=sig.u, col="blue")


## ------------------------------------------------------------------------
marginal = res$marginals.hyperpar$`Rho for t`

plot(marginal, type="l", xlab = expression(sigma), ylab = "Probability density",
     main = "Correlation of AR1 component", xlim=c(0.87, 1))
    
xvals = seq(0.85, 1, length.out=1000)
lines(xvals, 5*inla.pc.dcor1(xvals, hyper.ar1$theta2$param[1],
                              hyper.ar1$theta1$param[2]), lty='dashed')
  
abline(v=rho, col="blue")


## ------------------------------------------------------------------------
plot(res$summary.random$year$mean)
points(seas, col="blue")

