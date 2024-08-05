## ----setup, include=FALSE--------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
rm(list=ls())


## ---- warning=FALSE, message=FALSE-----------------------------------------------
library(INLA)
library(fields)
library(viridisLite)
library(gtools)
library(mapdata)
## Load functions
source("functions-rgeneric-121-march2020.R")


## --------------------------------------------------------------------------------
## Number of time points used
## Subsets the data
t.max = 4
## Max edge in spatial mesh
## Small number makes algorithm very slow
## ex: 2, 1.5
max.e = 2


## --------------------------------------------------------------------------------
## Prior model information
M = list()
for (i in 1:3) M[[i]] = list()
## Posteriors
fits = list()
## stacks (see later)
stack = list()


## --------------------------------------------------------------------------------
load("data/btopic133/tmed-eur.Rdata")


## --------------------------------------------------------------------------------
str(df.eur)


## --------------------------------------------------------------------------------
df1 = data.frame(locx = df.eur$longitude, locy = df.eur$latitude, 
                 time = df.eur$week, y = df.eur$tmed)


## --------------------------------------------------------------------------------
plot(df1$locx, df1$locy, col="blue",
     xlab="longitude", ylab="latitude", main="Observation locations")
map("worldHires", add=TRUE, col = grey(.5))


## --------------------------------------------------------------------------------
hist(df1$y, breaks=100)


## --------------------------------------------------------------------------------
## Transform with our chosen function
alt = log(df.eur$altitude+10)
## We group this into classes
alt = inla.group(alt, n=200)
## We add this to the dataframe
## We note that it has been transformed in some way by writing .t
df1$altitude.t = alt


## --------------------------------------------------------------------------------
hist(df1$altitude.t, breaks=100)


## --------------------------------------------------------------------------------
if (T) { 
  stopifnot(t.max <= 20)
  df2 = df1[df1$time <= t.max, ]
} else {
  ## Alternative: use every 3rd week and recode time variable
  ## This results in weaker time dependency
  stopifnot(t.max < 7)
  df2 = df1[df1$time %in% ((1:t.max)*3), ]
  df2$time = df2$time/3
}

summary(df2)


## --------------------------------------------------------------------------------
plot(df2$time, df2$y)


## --------------------------------------------------------------------------------
for (i in 1:t.max) {
  df.temp = df2[df2$t == i, ]
  print(mean(df.temp$y))
}


## --------------------------------------------------------------------------------
plot(df2$altitude.t, df2$y)


## --------------------------------------------------------------------------------
for (i in 1:3) {
  df.temp = df2[df2$time == i, ]
  quilt.plot(x = df.temp$locx, y = df.temp$locy, z = df.temp$y,
             nx = 40, ny = 40, col = plasma(101),
             main = paste("Plot data in year", i), zlim=range(df2$y))
  
}


## --------------------------------------------------------------------------------
plot(df2$locy, df2$y)
abline(lm(df2$y ~ df2$locy), col="blue")


## --------------------------------------------------------------------------------
hyper.iid.fix = list(prec = list(initial = -2*log(0.1), fixed=T))
hyper2 = list(prec = list(prior="pc.prec", param=c(0.5,0.5)))
hyper.ar1.rho = list(rho = list(prior = "pc.cor1", param=c(0.95,0.5)))


## --------------------------------------------------------------------------------
## Formula
form.lmm = y ~  locy + f(time, model = "iid", hyper = hyper2) + f(altitude.t, model = "rw1", hyper = hyper.iid.fix, scale.model = T)

## Fit model
fits[[1]] = inla(form.lmm,
                 family="gaussian",
                 data=df2,
                 num.threads = 3,
                 control.predictor = list(compute=T),
                 control.inla = list(int.strategy = "eb"))


## --------------------------------------------------------------------------------
summary(fits[[1]])


## --------------------------------------------------------------------------------
resid = df2$y - fits[[1]]$summary.linear.predictor$mean
hist(resid, breaks = 200)


## --------------------------------------------------------------------------------
for (i in 1:3) {
  df.temp = df2[df2$t == i, ]
  resid.temp = resid[df2$t == i]
  quilt.plot(x = df.temp$locx, y = df.temp$locy, z = resid.temp, 
             nx = 40, ny = 40, col = plasma(101),
             main = paste("Residuals year", i), zlim = range(resid))
  
}


## --------------------------------------------------------------------------------
mesh.t = inla.mesh.1d(1:t.max)


## --------------------------------------------------------------------------------
mesh.s = inla.mesh.2d(loc = cbind(df2$locx, df2$locy), 
                      max.edge=max.e*c(1, 2),
                      cutoff=max.e/5)


## --------------------------------------------------------------------------------
plot(mesh.s)
axis(1)
axis(2)


## --------------------------------------------------------------------------------
mesh.s$n*mesh.t$n


## --------------------------------------------------------------------------------
## Model component in space
mco.space = inla.spde2.pcmatern(mesh = mesh.s, 
                                prior.range = c(100, .5), prior.sigma = c(1, .5))


## --------------------------------------------------------------------------------
A.s1 = inla.spde.make.A(mesh = mesh.s, loc = cbind(df2$locx, df2$locy))


## --------------------------------------------------------------------------------
iset = inla.spde.make.index('i', n.spde = mesh.s$n, n.group = t.max)
A.st = inla.spde.make.A(mesh = mesh.s, loc = cbind(df2$locx, df2$locy), group = df2$time) 


## --------------------------------------------------------------------------------
stack[[2]] = inla.stack(
  data = list(y = df2$y), 
  effects = list(iset, m = rep(1, nrow(df2)), 
                 time = df2$time, altitude.t = df2$altitude.t,
                 s = 1:mesh.s$n),
  A = list(A.st, 1, 1, 1, A.s1), 
  tag = 'stdata') 


## --------------------------------------------------------------------------------
M[[2]]$shortname = "Separable"


## --------------------------------------------------------------------------------
## The separable model
hyper22 = list(prec = list(initial = -2*log(5), fixed=T))

## The formula
## It is possible to add a spatial random effect (uncomment f(s))
form.sep = y ~ -1 + m  + f(i, model = mco.space, group = i.group, 
                     control.group = list(model="ar1", hyper = hyper.ar1.rho)) + f(time, model = "iid", hyper=hyper22) + f(altitude.t, model = "rw1", hyper=hyper.iid.fix, scale.model = T) #+ f(s, model=mco.space)

M[[2]]$formula = form.sep


## --------------------------------------------------------------------------------
## Rgeneric object containing needed variables
## Mesh in space and time
## Lambdas for exponential prior on transformed hyper-param (1/rt, 1/rs and sig)
rgen.obj = list(mesh.space = mesh.s,
                mesh.time = mesh.t,
                lambdas = c(1,1,5))

## Nonsep model definition
nm = mesh.s$n*mesh.t$n

## The non-separable random effect / random field
mco.nonsep = inla.rgeneric.define(
  model = stmodel121.interpret, debug = FALSE, n = nm, obj = rgen.obj)


## --------------------------------------------------------------------------------
## Special: Since the index is very structured we can do this
## We know that iset is in the right order 
## (year 1 then year 2 etc, with the whole spatial mesh each time)
i.nonsep = 1:(mesh.s$n * mesh.t$n)

stack[[3]] = inla.stack(
  data = list(y = df2$y), 
  effects = list(i.nonsep = i.nonsep, m = rep(1, nrow(df2)), 
                 time = df2$time, altitude.t = df2$altitude.t,
                 s = 1:mesh.s$n),
  A = list(A.st, 1, 1, 1, A.s1), 
  tag = 'stdata') 


## --------------------------------------------------------------------------------
M[[2]]$shortname = "Non-separable"


## --------------------------------------------------------------------------------
mco.nonsep.fix = inla.rgeneric.define(model = stmodel121.interpret, 
                                      debug = FALSE, n = nm, obj = rgen.obj)

## The formula
## It is possible to add a spatial random effect (uncomment the f(s))
form.nonsep = y ~ -1 + m + f(i.nonsep, model = mco.nonsep.fix, n = nm)   + f(time, model = "iid", hyper = hyper22) + f(altitude.t, model = "rw1", hyper = hyper.iid.fix, scale.model = T) #+ f(s, model=mco.space)

M[[3]]$formula = form.nonsep


## --------------------------------------------------------------------------------
M[[2]]$init = c(0.404 , 2.207 , 1.813 , 3.804)
M[[3]]$init = c(0.435 , 8.012 , 2.944 , 2.553)


## --------------------------------------------------------------------------------
## Fit model 2 and 3
for (i in 2:3){
  print(paste("Running:  ", i))
  stk = stack[[i]]
  fits[[i]] = inla(M[[i]]$formula,
                   family="gaussian",
                   data=inla.stack.data(stk),
                   control.predictor=list(A=inla.stack.A(stk), compute=T),
                   #verbose=T,
                   num.threads = 3,
                   control.inla = list(int.strategy = "eb"),
                   control.mode = list(restart=T, theta=M[[i]]$init),
                   control.compute = list(dic=TRUE, cpo=TRUE, waic=T, 
                           mlik=T, return.marginals=F, config=T, 
                           openmp.strategy="default", smtp="taucs")
)
}


## --------------------------------------------------------------------------------
for (i in 1:length(M)){
    print(round(fits[[i]]$cpu.used[4],2))
}


## --------------------------------------------------------------------------------
for (i in 1:length(M)){
    print(paste(round(fits[[i]]$internal.summary.hyperpar$mean, 3), collapse = " , "))
}


## --------------------------------------------------------------------------------
## Comparison of the Precision for Gaussian observations
fits[[2]]$summary.hyperpar
fits[[3]]$summary.hyperpar

## To compare the 3 other hyper-parameters, we need to transform the
## hyperparameters of the non-separable model as follows
data.frame(var=c("Range.s", "Stdev", "Range.t"), exp(fits[[3]]$summary.hyperpar[c(3,4,2), c(4,3,5)]))


## --------------------------------------------------------------------------------
#summary(fits[[2]])
#summary(fits[[3]])


## --------------------------------------------------------------------------------
plot(df2$time, df2$y)
points(fits[[2]]$summary.random$time$mean, col="blue", pch=19, cex=2)


## --------------------------------------------------------------------------------
plot(fits[[2]]$summary.random$altitude$mean)


## --------------------------------------------------------------------------------
## Compare fitted values and observations.
# plot(fits[[2]]$summary.linear.predictor$mean[1:nrow(df2)], df2$y)


## --------------------------------------------------------------------------------
rf.st.sep = fits[[2]]$summary.random$i$mean + fits[[2]]$summary.fixed$mean[1]
rf.st.nonsep = fits[[3]]$summary.random$i$mean + fits[[3]]$summary.fixed$mean[1]


## --------------------------------------------------------------------------------
## Define default plotting window
local.xlim = range(df2$locx)
local.ylim = range(df2$locy)

local.plot.field = function(field, mesh=mesh.s, 
                            timepoint=1, xlim, ylim, ...){
  ## Here: mesh is spatial mesh only
  # Possible error when using a wrong spatial mesh
  stopifnot((length(field) %% mesh$n) == 0 )
  
  field = field[1:mesh$n + (timepoint-1)*mesh$n]
  # - only use the relevant part of the incoming vector
  
  # Choose plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim = local.xlim
  if (missing(ylim)) ylim = local.ylim
  
  # Project the mesh onto a 300x300 grid
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  
  # Do the projection 
  field.proj = inla.mesh.project(proj, field)
  
  # Plot it
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), col = plasma(64),
             xlim = xlim, ylim = ylim, ...)  
}


## --------------------------------------------------------------------------------
if(!is.null(fits[[2]]$summary.random$s)) {
  local.plot.field(fits[[2]]$summary.random$s$mean)
}


## --------------------------------------------------------------------------------
if(!is.null(fits[[3]]$summary.random$s)) {
  local.plot.field(fits[[3]]$summary.random$s$mean)
}


## --------------------------------------------------------------------------------
## Change this to any year
time.k = 3


## --------------------------------------------------------------------------------
local.plot.field(rf.st.sep, timepoint = time.k, 
                 main = paste("Separable model at time =", i))
points(df2$locx, df2$locy, col=adjustcolor("black", alpha.f = 0.01))
map("worldHires", add=TRUE, col = grey(.5))


## --------------------------------------------------------------------------------
local.plot.field(rf.st.nonsep, timepoint = time.k, 
                 main = paste("Nonsep model at time =", time.k))
points(df2$locx, df2$locy, col=adjustcolor("black", alpha.f = 0.01))
map("worldHires", add=TRUE, col = grey(.5))


## --------------------------------------------------------------------------------
local.plot.field(rf.st.nonsep-rf.st.sep, timepoint = time.k)
points(df2$locx, df2$locy, col=adjustcolor("black", alpha.f = 0.01))
map("worldHires", add=TRUE, col = grey(.5))


## --------------------------------------------------------------------------------
zlim.temp = range(rf.st.nonsep) #c(-9, 16)
for (i in 1:3) {
  local.plot.field(rf.st.nonsep, timepoint = i, zlim=zlim.temp,
                   main=paste("Nonsep model at time =", i))
  map("worldHires", add=TRUE, col = grey(.5))
}

