## ----setup, include=FALSE--------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- include=FALSE--------------------------------------------------------------
rm(list=ls())
options(width=70, digits=2)


## ---- warning=FALSE, message=FALSE-----------------------------------------------
library(INLA)


## --------------------------------------------------------------------------------
data(Seeds)


## --------------------------------------------------------------------------------
# Run: ?Seeds
head(Seeds)
# - r is the number of seed germinated (successes)
# - n is the number of seeds attempted (trials)
# - x1 is the type of seed
# - x2 is the type of root extract
# - plate is the numbering of the plates/experiments


## --------------------------------------------------------------------------------
df = data.frame(y = Seeds$r, Ntrials = Seeds$n, Seeds[, 3:5])


## --------------------------------------------------------------------------------
summary(df)
table(df$x1)
table(df$x2)


## --------------------------------------------------------------------------------
plot(df)


## --------------------------------------------------------------------------------
family1 = "binomial"
control.family1 = list(control.link=list(model="logit"))
# number of trials is df$Ntrials


## --------------------------------------------------------------------------------
hyper1 = list(theta = list(prior="pc.prec", param=c(1,0.01)))
formula1 = y ~ x1 + x2 + f(plate, model="iid", hyper=hyper1)


## --------------------------------------------------------------------------------
res1 = inla(formula=formula1, data=df, 
            family=family1, Ntrials=Ntrials, 
            control.family=control.family1)


## --------------------------------------------------------------------------------
summary(res1)


## --------------------------------------------------------------------------------
# Run: str(res1, 1)


## --------------------------------------------------------------------------------
res1$summary.random$plate


## --------------------------------------------------------------------------------
plot(1:nrow(df), res1$summary.random$plate$`0.97`, col="red", ylim=c(-1,1),
     xlab="measurement number", ylab = "quantiles")
points(1:nrow(df), res1$summary.random$plate$`0.5quant`)
points(1:nrow(df), res1$summary.random$plate$`0.02`, col="blue")


## --------------------------------------------------------------------------------
m.beta1 = inla.tmarginal(fun = function(x) x, marginal = 
                           res1$marginals.fixed$x1)
# - this transformation is the identity (does nothing)
# - m.beta1 is the marginal for the coefficient in front of the x1 covariate
plot(m.beta1, type="l", xlab = expression(beta[1]), ylab = "Probability density")


## --------------------------------------------------------------------------------
m.sigma = inla.tmarginal(fun = function(x) exp(-1/2*x), marginal = 
                           res1$internal.marginals.hyperpar$`Log precision for plate`)
# - m.sigma is the marginal for the standard deviation parameter in the iid random effect
plot(m.sigma, type="l", xlab = expression(sigma[iid]), ylab = "Probability density")


## --------------------------------------------------------------------------------
m.plate.7 = inla.tmarginal(fun = function(x) x, marginal = 
                           res1$marginals.random$plate$index.7)

# - m.plate.7 is one of the marginals for the parameters beloning to the plate iid effect
# - it is number 7, which corresponds to plate=7, which is our 7th row of data
plot(m.plate.7, type="l", xlab = "marginal plate nr 7", ylab = "Probability density")


## --------------------------------------------------------------------------------
plot(density(res1$summary.random$plate$mean))
lines(0+c(-2, 2)*res1$summary.hyperpar$`0.5quant`^(-0.5) , c(0,0), col="blue")
# - draw a blue line for plus/minus 2 sigma


## --------------------------------------------------------------------------------
df2 = rbind(df, c(NA, 1, 0, 0, 22))
tail(df2)


## --------------------------------------------------------------------------------
res.pred = inla(formula=formula1, data=df2, 
            family=family1, Ntrials=Ntrials, 
            control.predictor = list(compute=T, link = 1),
            # - to get the posterior of the predictor and fitted values
            control.family=control.family1)


## --------------------------------------------------------------------------------
res.pred$summary.fitted.values[22, ]
# - this is the inv.logit(eta_i), namely p_i


## --------------------------------------------------------------------------------
control.family2 = list(control.link=list(model="probit"))
res1 = inla(formula=formula1, data=df, 
            family=family1, Ntrials=Ntrials, 
            control.predictor = list(compute=T, link=1),
            # - to get the posterior of the predictor and fitted values
            control.family=control.family1)
res2 = inla(formula=formula1, data=df, 
            family=family1, Ntrials=Ntrials, 
            control.predictor = list(compute=T, link=1),
            # - to get the posterior of the predictor and fitted values
            control.family=control.family2)


## --------------------------------------------------------------------------------
a = data.frame(y=df$y, ps=df$y/df$Ntrials,
               r1eta = res1$summary.linear.predictor$mean,
           r2eta = res2$summary.linear.predictor$mean,
           r1fit = res1$summary.fitted.values$mean,
           r2fit = res2$summary.fitted.values$mean)
round(a, 2)

