## ----setup, include=FALSE--------------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## --------------------------------------------------------------------------------
library(ggplot2)
set.seed(2021)


## --------------------------------------------------------------------------------
b = 0.0076
lambda = 4.75
tau = seq(0.0001, 200, length.out=1E4)
pri.tau.old = b*exp(-b*tau)
pri.tau.new = lambda/2 * tau^(-3/2) *exp(-lambda * tau^(-1/2))
plot(tau, pri.tau.old, type="l", ylim=c(0, max(c(pri.tau.old, pri.tau.new))), ylab="density")
lines(tau, pri.tau.new, col="blue")


## --------------------------------------------------------------------------------
sigma = seq(0.0001, 0.6, length.out=1E4)
pri.sig.old = 2*b*sigma^(-3) *exp(-b*sigma^(-2))
pri.sig.new = lambda*exp(-lambda*sigma)
plot(sigma, pri.sig.old, type="l", ylim=c(0, max(c(pri.sig.old, pri.sig.new))), ylab="density")
lines(sigma, pri.sig.new, col="blue")


## ---- message=FALSE--------------------------------------------------------------
library(INLA)


## --------------------------------------------------------------------------------
corvals = seq(0.01, 0.99, length.out = 1000)
d = inla.pc.dcor1(corvals, lambda=1)
plot(corvals, d, type="l")

