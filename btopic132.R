## ----setup, include=FALSE-----------------------------------------
knitr::opts_chunk$set(echo = TRUE)
rm(list=ls())


## ---- warning=FALSE, message=FALSE--------------------------------
library(INLA)
library(fields)
library(viridisLite)
library(gtools)
## Load functions
source("functions-rgeneric-121-march2020.R")
## Seed
set.seed(20200119)


## -----------------------------------------------------------------
## Range in time
## For nonseparable 
## (we use a factor of 1.8 modification when comparing to separable, see paper)
range.t = 3.5

## Range in space
range.s = 5

## Max edge in spatial mesh
## Small numbers makes algorithm very slow
## Should be less than 1/4 of range.s, preferably much smaller
## ex: 1.5, 1, 0.7, 0.5
max.e = 1

## Number of timepoints used
## Must be 2 or greater
t.max = 5


## -----------------------------------------------------------------
mesh.t = inla.mesh.1d(1:t.max)


## -----------------------------------------------------------------
fake.locations = matrix(c(0,0,10,10, 0, 10, 10, 0), nrow = 4, byrow = T)
mesh.s = inla.mesh.2d(loc = fake.locations, max.edge=max.e*c(1, 2))


## -----------------------------------------------------------------
plot(mesh.s)


## -----------------------------------------------------------------
## Approx rho for decay in 6 time points
approx.rho = 0.13^(1/6)
theta2 = logit(approx.rho)
N = t.max 
## old syntax: Q = INLA:::inla.extract.Q('i', formula = y ~ f(i, model='ar1', hyper = list(prec=list(initial=0, fixed=T), theta2=list(initial=theta2, fixed=T))), data = data.frame(y=1:N, i=1:N))
res1 = inla(formula = y ~ f(i, model='ar1', hyper = list(prec=list(initial=0, fixed=T), theta2=list(initial=theta2, fixed=T))), data = data.frame(y=1:N, i=1:N), control.compute=list(config=T))
## Does not work: Q = INLA:::inla.extract.Q("i", result = res1)
if (T) {
  ## New way, but bad code:
  ## Lowest row of matrix a 0
  res1$misc$configs$contents$tag
  result=res1
  what = "i"
  conf <- result$misc$configs
  k=6
  Q1 <- conf$config[[k]]$Qprior
  Q = t(Q1)+Q1
  diag(Q) = diag(Q)/2
  Q = Q[1:5, 1:5]
  
}
Q
solve(Q)


## -----------------------------------------------------------------
mco.space = inla.spde2.pcmatern(mesh = mesh.s, prior.range = c(5, .5), prior.sigma = c(1, .5))


## -----------------------------------------------------------------
Qsep.s = inla.spde2.precision(spde = mco.space, theta = log(c(range.s,1)))


## -----------------------------------------------------------------
## Gaussian noise
sig.eps = 0.01
## Seed used in C code
inla.seed = sample(1E12, 1)
## Sample with INLA
u.sim = inla.qsample(n=1, Qsep.s, seed = inla.seed, num.threads=1)[, 1]
u.sim = u.sim - mean(u.sim)
sim.noise = rnorm(length(u.sim), 0, 1) * sig.eps
## st is spacetime index
df1 = data.frame(y=u.sim+sim.noise,
                u.sim = u.sim, sim.noise=sim.noise,
                year=1)
summary(df1)


## -----------------------------------------------------------------
## df2: Augment df1 with the needed prediction locations
temp.na = rep(NA, (t.max-1)*mesh.s$n)
df2 = data.frame(y=c(df1$y, temp.na), 
                 locx = rep(mesh.s$loc[, 1], t.max),
                 locy = rep(mesh.s$loc[, 2], t.max),
                 year = rep(1:t.max, each=mesh.s$n)
)

## Our final dataframe
summary(df2)


## -----------------------------------------------------------------
## Rgeneric object containing needed variables
## Mesh in space and time
## Lambdas for exponential prior on ransformed hyper-param (1/rt, 1/rs and sig)
rgen.obj = list(mesh.space = mesh.s,
                mesh.time = mesh.t,
                lambdas = c(1,1,1))

## Nonsep model definition
nm = mesh.s$n*mesh.t$n

## The non-separable random effect / random field
## We use the function loaded in the beginning of the document
mco.nonsep = inla.rgeneric.define(
  model=stmodel121.interpret, debug=FALSE, n=nm, obj=rgen.obj)


## -----------------------------------------------------------------
local.plot.field = function(field, mesh=mesh.s, 
                            timepoint=1, xlim, ylim, ...){
  ## Here: mesh is spatial mesh only
  # Possible error when using a wrong spatial mesh
  stopifnot((length(field) %% mesh$n) ==0 )
  
  field = field[1:mesh$n + (timepoint-1)*mesh$n]
  # - only use the relevant part of the incoming vector
  
  # Choose plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim = c(0, 10)
  if (missing(ylim)) ylim = c(0, 10)
  
  # Project the mesh onto a 300x300 grid
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  
  # Do the projection 
  field.proj = inla.mesh.project(proj, field)
  
  # Plot it
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), col = plasma(64),
             xlim = xlim, ylim = ylim, ...)  
}


## -----------------------------------------------------------------
local.plot.field(field = df1$u.sim, mesh = mesh.s)


## -----------------------------------------------------------------
M = list()
for (i in 1:2) M[[i]] = list()


## -----------------------------------------------------------------
iset = inla.spde.make.index('i', n.spde = mesh.s$n, n.group = t.max)
A = inla.spde.make.A(mesh = mesh.s, loc = cbind(df2$locx, df2$locy), group = df2$year) 
sum(A-Diagonal(nrow(A)))


## -----------------------------------------------------------------
stack = list()
stack[[1]] = inla.stack(
  data = list(y = df2$y), 
  A = list(A, 1), 
  effects = list(iset, m = rep(1, nrow(df2))),
  tag = 'stdata') 


## -----------------------------------------------------------------
## Special: Since the index is very structured we can use the same A matrix here
## We know that iset is in the right order 
## (year 1 then year 2 etc, with the whole spatial mesh each time)
i.nonsep = 1:nm

stack[[2]] = inla.stack(
  data = list(y = df2$y), 
  A = list(A, 1), 
  effects = list(i.nonsep=i.nonsep, m = rep(1, nrow(df2))),
  tag = 'stdata') 


## -----------------------------------------------------------------
## We want to fix the autocorrelation in time
## The theta2 as defined before is the INLA internal parametrisation of rho
hyper.ar1.rho = list(rho = list(initial=theta2, fixed=TRUE))
form.sep = y ~ -1 + m + f(i, model = mco.space, group = i.group, 
                     control.group = list(model="ar1", hyper = hyper.ar1.rho))
M[[1]]$formula = form.sep


## -----------------------------------------------------------------
## We need to fix the temporal range in the non-separable model
rgen.obj2 = rgen.obj
rgen.obj2$fixed.theta = c(log(range.t), NA, NA)
mco.nonsep.fix = inla.rgeneric.define(model=stmodel121.interpret, debug=FALSE, n=nm, obj=rgen.obj2)

form.nonsep = y ~ -1 + m + f(i.nonsep, model=mco.nonsep.fix, n=nm)
M[[2]]$formula = form.nonsep


## -----------------------------------------------------------------
print(M)


## -----------------------------------------------------------------
## WARNING: Set these variables to NULL if you change the model in any way!
M[[1]]$init = c(1.775 , 0.153)
M[[2]]$init = c(1.11 , -0.056)

fits = list()


## -----------------------------------------------------------------
for (i in 1:length(M)){
  print(paste("Running:  ", i))
  stk = stack[[i]]
  fits[[i]] = inla(M[[i]]$formula,
                   family="gaussian",
                   control.family = list(hyper = list(prec = list(
                     initial = log(sig.eps)*-2, fixed=T)
                   )),
                   data=inla.stack.data(stk),
                   control.predictor=list(A=inla.stack.A(stk), compute=T),
                   verbose=F,
                   num.threads = 3,
                   control.inla = list(int.strategy = "eb"),
                   control.mode = list(restart=T, theta=M[[i]]$init),
                   control.compute = list(dic=TRUE, cpo=TRUE, waic=T, 
                           mlik=T, return.marginals=F, config=T, 
                           openmp.strategy="default", smtp="taucs")
)
    print(round(fits[[i]]$cpu.used[4],2))
}


## -----------------------------------------------------------------
## Check what we can set the initial values to
for (i in 1:length(M)){
    print(paste(round(fits[[i]]$internal.summary.hyperpar$mean, 3), collapse = " , "))
}


## -----------------------------------------------------------------
## Comparison
fits[[1]]$summary.hyperpar[ ,c(4,3,5)]

## This only works in this specific case, when we fixed the first hyper-parameter
data.frame(var=c("Range", "Stdev"), exp(fits[[2]]$summary.hyperpar[c(1,2), c(4,3,5)]))


## -----------------------------------------------------------------
pred.sep = fits[[1]]$summary.random$i$mean + fits[[1]]$summary.fixed$mean[1]
pred.nonsep = fits[[2]]$summary.random$i$mean + fits[[2]]$summary.fixed$mean[1]


## -----------------------------------------------------------------
## zlim for plots
local.zlim = c(-2, 3)
## Yearly reduction in zlim
zlim.disc = approx.rho # discount factor


## -----------------------------------------------------------------
local.plot.field(pred.sep, timepoint = 1, zlim=local.zlim)


## -----------------------------------------------------------------
local.plot.field(pred.nonsep, timepoint = 1, zlim=local.zlim)


## -----------------------------------------------------------------
## Difference between predictions.
local.plot.field(pred.nonsep-pred.sep, timepoint = 1, zlim=c(-1,1))


## -----------------------------------------------------------------
local.plot.field(pred.sep, timepoint = 2, zlim=local.zlim*zlim.disc)


## -----------------------------------------------------------------
local.plot.field(pred.nonsep, timepoint = 2, zlim=local.zlim*zlim.disc)


## -----------------------------------------------------------------
local.plot.field(pred.nonsep-pred.sep, timepoint = 2, zlim=c(-1,1))


## -----------------------------------------------------------------
if(t.max>=5) local.plot.field(pred.sep, timepoint = 5, zlim=local.zlim*zlim.disc^5)


## -----------------------------------------------------------------
if(t.max>=5) local.plot.field(pred.nonsep, timepoint = 5, zlim=local.zlim*zlim.disc^5)


## -----------------------------------------------------------------
if(t.max>=5) local.plot.field(pred.nonsep-pred.sep, timepoint = 5, zlim=c(-1,1))

