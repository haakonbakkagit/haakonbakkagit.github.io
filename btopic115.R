## ----setup, include=FALSE-------------------------------------------
rm(list = ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE----------------------------------
library(INLA)
library(shiny)


## -------------------------------------------------------------------
temp = read.csv("data/harmonised-unemployment-rates-mo.csv")
n = nrow(temp)
data = data.frame(y = temp[,2], t=1:n)
dates <- temp[,1]


## -------------------------------------------------------------------
plot(dates, data$y, lwd=2,
    xlab='month', ylab='Unemployment Rates')
lines(dates,data$y)
abline(h=2*(-8:9), lty=2, col=gray(.5))


## -------------------------------------------------------------------
family <- "gaussian"


## -------------------------------------------------------------------
formula1 <- y~ f(t,model='ar1')


## -------------------------------------------------------------------
res1 <- inla(formula=formula1,data=data,family=family)


## -------------------------------------------------------------------
summary(res1)


## -------------------------------------------------------------------
plot(res1$summary.random$t[ ,"mean"]+res1$summary.fixed$mean[1],ylab="fitting result",type="l")
points(data$y, col="blue")


## -------------------------------------------------------------------
hyper2 = list(theta1=list(initial=0.5, fixed=T))
formula2 <- y~ f(t,model='ar1',hyper=hyper2)
  
res2 <- inla(formula=formula2,data=data,family=family)
plot(data$y, col="blue",
     ylab="fitting result")
lines(res2$summary.random$t[ ,"mean"]+res2$summary.fixed$mean[1])


## -------------------------------------------------------------------
family <- "gaussian"

hyper3 <- list(theta1 = list(prior="pc.prec", param=c(0.06, 0.008)),
                    theta2 = list(prior="pc.cor1", param=c(0.9, 0.9)) )
formula3 <- y~ f(t,model='ar1',hyper=hyper3)
res3 <- inla(formula=formula3,data=data,family=family,
             control.predictor = list(compute=T))

plot(data$y, col="blue",
     ylab="fitting result")
lines(res3$summary.random$t[ ,"mean"]+res3$summary.fixed$mean[1])


## -------------------------------------------------------------------
plot(1:n, res3$summary.random$t$`0.97`, col="red", type="l",
   ylim=c(-6,6),xlab="measurement number", ylab = "quantiles")
lines(1:n, res3$summary.random$t$`0.5quant`)
lines(1:n, res3$summary.random$t$`0.02`, col="blue")


## -------------------------------------------------------------------
m.sigma = inla.tmarginal(fun=function(x)x^(-0.5),marginal = 
                           res3$marginals.hyperpar$`Precision for t`)
plot(m.sigma, type="l", xlab = expression(sigma), ylab = "Probability density")
xvals = seq(0.5, 1.5, length.out=1000)
lambda=-log(0.008)/0.06
lines(xvals, 1e30*lambda*exp(-lambda*xvals), lty='dashed')


## -------------------------------------------------------------------
m.rho <- inla.tmarginal(fun=function(x)x,marginal = 
                            res3$marginals.hyperpar$`Rho for t`)
plot(m.rho, type="l", xlab = expression(rho), ylab = "Probability density")
xvals = seq(0.5, 1, length.out=1000)
lines(xvals, 5*inla.pc.dcor1(xvals, 0.9 , 0.9 , log = FALSE), lty='dashed')


## -------------------------------------------------------------------
m.t.70 <- inla.tmarginal(fun = function(x) x, marginal = 
                             res3$marginals.random$t$index.70)
# - m.t.70 is one of the marginals for the parameters beloning to the plate iid effect
# - it is number 70, which corresponds to plate=70, which is our 70th row of data
plot(m.t.70, type="l", xlab = "marginal t nr 70", ylab = "Probability density")


## -------------------------------------------------------------------
m.theta3 <- inla.tmarginal(fun=function(x)x,marginal = 
                            res3$marginals.fixed$`(Intercept)`)
plot(m.theta3, type="l", xlab = "intercept", ylab = "Probability density")


## -------------------------------------------------------------------
hist(data$y - res3$summary.fitted.values$mean, breaks = 50)

