## ----setup, include=FALSE-------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE----------------------------------
library(INLA)


## -------------------------------------------------------------------
temp = read.csv("data/harmonised-unemployment-rates-mo.csv")
n = nrow(temp)
data = data.frame(y = temp[,2], t=1:n)
dates <- temp[,1]

df = data.frame(data, dates)

summary(df)


## -------------------------------------------------------------------
plot(df$dates, df$y, lwd=2, col="blue", xlab='Month', ylab='Percentage',
     main = "Unemployment Norwegian Females (standardised)")
lines(dates,data$y)
abline(h=2*(-8:9), lty=2, col=gray(.5))


## -------------------------------------------------------------------
formula = y~ f(t,model='ar1')


## -------------------------------------------------------------------
res = inla(formula=formula,data=df,family="gaussian",
           control.predictor=list(compute=TRUE))


## -------------------------------------------------------------------
summary(res)


## -------------------------------------------------------------------
str(res$summary.random$t)


## -------------------------------------------------------------------
str(res$marginals.fixed)


## -------------------------------------------------------------------
str(res$marginals.hyperpar, 1)


## -------------------------------------------------------------------
plot(df$y, col="blue", main="Fitting result 1", xlab=NA, ylab=NA)
lines(res$summary.fitted.values$mean)
lines(res$summary.fitted.values$`0.02`,col="grey")
lines(res$summary.fitted.values$`0.97`,col="grey")


## ---- eval=FALSE----------------------------------------------------
## inla.doc("pc.prec")
## inla.doc("pc.cor1")


## -------------------------------------------------------------------
hyper.ar1 = list(theta1 = list(prior="pc.prec", param=c(0.02, 0.5)),
                 theta2 = list(prior="pc.cor1", param=c(0.9, 0.5)))


## -------------------------------------------------------------------
hyper.family = list(theta = list(prior="pc.prec", param=c(3, 0.5)))


## -------------------------------------------------------------------
formula2 <- y~ f(t,model='ar1', hyper=hyper.ar1)

res2 <- inla(formula=formula2,data=df,family="gaussian",
             control.predictor=list(compute=TRUE),
             control.family = list(hyper = hyper.family))


## -------------------------------------------------------------------
plot(df$y, col="blue", main="Fitting result 2", xlab=NA, ylab=NA)
lines(res2$summary.fitted.values$mean)
lines(res2$summary.fitted.values$`0.02`,col="grey")
lines(res2$summary.fitted.values$`0.97`,col="grey")

