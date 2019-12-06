rm(list = ls())
library(INLA); library(sp); library(fields)
source('functions-barrier-DT.R')


### DESCRIPTION ###
# This simulates the first troll dataset
# THIS FILE IS NOT USED YET
#  Instead I use the Ca20 data

### INPUT ###
range = 2
N = 100
set.seed(2017)
set.inla.seed = 2017


### MESH ###
poly = bakka.square.polygon(xlim = c(0, 10), ylim=c(0, 10), ret.SP=T)

max.edge = range/10
bound.outer = 3
mesh = inla.mesh.2d(boundary = poly,
                    max.edge = c(1,5)*max.edge,
                    cutoff = max.edge/2,
                    offset = c(max.edge, bound.outer))

mesh$n
plot(mesh)

### PRECISION ###
spde = inla.spde2.pcmatern(mesh, prior.range = c(5, .5), prior.sigma = c(.5, .5))
# - ignore the priors, they are not used at all (in this topic)
Q = inla.spde2.precision(spde, theta = c(log(range),0))
u = inla.qsample(n=1L, Q=Q, seed = set.inla.seed)

### Random locations ###
locx = runif(n = N, min = 0, max=10)
locy = runif(n = N, min = 0, max=10)
A = inla.spde.make.A(mesh = mesh, loc = cbind(locx, locy))

plot(locx, locy)

Au = A %*% u

### Temperature ###
beta = function(x) ((1-(x/10)^2))

temp = 10 + 20*cos(locx)

beta.temp = beta(temp)

### DATA ###
Au = as.numeric(Au)
y.n = beta.temp + Au + rnorm(N, mean = 0.1) + 10
#E = 10
#y = rpois(n = N, lambda = exp(E*y.n))
df = data.frame(y=y.n, locx, locy, temp)
summary(df)


