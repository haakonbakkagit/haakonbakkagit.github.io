## ----setup, include=FALSE-------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE----------------------------------
library(INLA); library(sp); library(fields)
library(geoR)
library(viridisLite)
# - for better colours
rm(list=ls())
options(width=70, digits=2)


## -------------------------------------------------------------------
data('ca20')
class(ca20)
summary(ca20)


## -------------------------------------------------------------------
df = data.frame(y = ca20$data, locx = ca20[[1]][ , 1], locy = ca20[[1]][ , 2], ca20[[3]])
spatial.scaling = 100
df$locx = (df$locx - min(df$locx))/spatial.scaling
df$locy = (df$locy - min(df$locy))/spatial.scaling
df$altitude = df$altitude - mean(df$altitude)
df$y = df$y-50


## -------------------------------------------------------------------
df$area1 = (df$area==1)*1
df$area2 = (df$area==2)*1


## -------------------------------------------------------------------
summary(df)

head(df)

cor(cbind(df[, 1:4], as.numeric(df[ , 5])))


## -------------------------------------------------------------------
plot(df$altitude, df$y)
abline(lm(df$y~df$altitude), col="red")


## -------------------------------------------------------------------
formula = list()
formula[[1]] = y ~ altitude
formula[[2]] = y ~ f(altitude, model="rw1", scale.model = T, hyper = list(prec = list(prior="pc.prec", param=c(1,0.01))))


## -------------------------------------------------------------------
prior.median.gaus.sd = 1
# - Think about this value
# - Remember sd(df$y)
family = 'gaussian'
control.family = list(hyper = list(prec = list(
  prior = "pc.prec", fixed = FALSE, param = c(prior.median.gaus.sd,0.5))))


## -------------------------------------------------------------------
res = list()
for (i in 1:2){
res[[i]] <- inla(formula[[i]], data=df,
            family = family,
            control.family = control.family,
            control.predictor = list(compute=T),
            control.inla = list(int.strategy='eb'),
            control.fixed = list(expand.factor.strategy='inla'))
}


## -------------------------------------------------------------------
summary(res[[1]])


## -------------------------------------------------------------------
res[[1]]$summary.fixed[, 1:5]


## -------------------------------------------------------------------
res[[2]]$summary.fixed[, 1:5]


## -------------------------------------------------------------------
plot(df$altitude, df$y)
effect = df$altitude*res[[1]]$summary.fixed$mean[2]
lines(df$altitude, effect, col="blue")


## -------------------------------------------------------------------
plot(df$altitude, df$y)
lines(sort(unique(df$altitude)), res[[2]]$summary.random$altitude$mean, col="blue")
lines(sort(unique(df$altitude)), res[[2]]$summary.random$altitude$`0.025quant`, col="blue")
lines(sort(unique(df$altitude)), res[[2]]$summary.random$altitude$`0.97quant`, col="blue")


## -------------------------------------------------------------------
plot(df$altitude, df$y)
effect = res[[1]]$summary.linear.predictor$mean
lines(df$altitude, effect, col="blue")


## -------------------------------------------------------------------
plot(df$altitude, df$y)
points(df$altitude, res[[2]]$summary.linear.predictor$mean, col="blue")
points(df$altitude, res[[2]]$summary.linear.predictor$`0.025quant`, col="blue")
points(df$altitude, res[[2]]$summary.linear.predictor$`0.97quant`, col="blue")

