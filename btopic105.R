## ----setup, include=FALSE--------------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE-----------------------------------------------
library(INLA) 
library(sp) 
library(fields)

set.seed(2016)
set.inla.seed = 2016


## --------------------------------------------------------------------------------
## Load data
load(file = "data/WebSiteData-Archipelago.RData")
# - if you have saved the file locally

## What is loaded
# - poly.water is our study area
# - df is our dataframe to be analysed
# - dat is the orginial dataframe
str(poly.water, 1)


## --------------------------------------------------------------------------------
max.edge = 0.95
mesh1 = inla.mesh.2d(loc=cbind(df$locx, df$locy),
                    max.edge = max.edge)
plot(mesh1, main="1st attempt"); points(df$locx, df$locy, col="blue")


## --------------------------------------------------------------------------------
max.edge = 0.95
# - as before
bound.outer = 4.6
# - as before
mesh4 = inla.mesh.2d(boundary = poly.water,
                     loc=cbind(df$locx, df$locy),
                    max.edge = c(1,5)*max.edge,
# - use 5 times max.edge in the outer extension/offset/boundary
                    cutoff = 0.06,
                    offset = c(max.edge, bound.outer))
plot(mesh4, main="4th attempt", lwd=0.5); points(df$locx, df$locy, col="red")


## --------------------------------------------------------------------------------
local.plot.field = function(field, mesh, xlim, ylim, ...){
  stopifnot(length(field) == mesh$n)
  # - error when using the wrong mesh
  if (missing(xlim)) xlim = poly.water@bbox[1, ] 
  if (missing(ylim)) ylim = poly.water@bbox[2, ]
  # - choose plotting region to be the same as the study area polygon
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 grid 
  #   for plots
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, ...)  
}


## --------------------------------------------------------------------------------
local.find.correlation = function(Q, location, mesh) {
  sd = sqrt(diag(inla.qinv(Q)))
  # - the marginal standard deviations
  
  A.tmp = inla.spde.make.A(mesh=mesh, loc = matrix(c(location[1], location[2]),1,2))
  # - create a fake A matrix, to extract the closest mesh node index
  id.node = which.max(A.tmp[1, ])
  # - index of the closest node
  
  print(paste('The location used was c(', 
              round(mesh$loc[id.node, 1], 4), ', ', 
              round(mesh$loc[id.node, 2], 4), ')' ))
  # - location of the closest node
  # - should be close to the location input
  # - sometimes used to plot a black dot
  
  ## Solve a matrix system to find the column of the covariance matrix
  Inode = rep(0, dim(Q)[1]); Inode[id.node] = 1
  covar.column = solve(Q, Inode)
  corr = drop(matrix(covar.column)) / (sd*sd[id.node])
  return(corr)
}


## ---- warning=FALSE--------------------------------------------------------------
spde = inla.spde2.pcmatern(mesh1, prior.range = c(5, .5), prior.sigma = c(.5, .5))
# - ignore the priors, they are not used at all (in this topic)
Q = inla.spde2.precision(spde, theta = c(log(4),log(1)))
# - theta: log range and log sigma (standard deviation parameter)
sd = sqrt(diag(inla.qinv(Q)))
local.plot.field(sd, mesh1)
points(df$locx, df$locy)


## --------------------------------------------------------------------------------
corr = local.find.correlation(Q, loc = c(16.4, 6.9), mesh1)
local.plot.field(corr, mesh1, zlim=c(0.1, 1))
points(16.52, 6.93)


## --------------------------------------------------------------------------------
corr = local.find.correlation(Q, loc = c(8, 7), mesh1)
local.plot.field(corr, mesh1, zlim=c(0.1, 1))
points(7.62, 6.77)


## ---- warning=FALSE--------------------------------------------------------------
spde = inla.spde2.pcmatern(mesh4, prior.range = c(5, .5), prior.sigma = c(.5, .5))
# - You can ignore the prior, we do not use that
Q = inla.spde2.precision(spde, theta = c(log(4),log(1)))
# - log range and log sigma (standard deviation)


## --------------------------------------------------------------------------------
sd = diag(inla.qinv(Q))
local.plot.field(sd, mesh4)
points(df$locx, df$locy)
plot(poly.water, add=T)


## --------------------------------------------------------------------------------
corr = local.find.correlation(Q, loc = c(7,10.3), mesh4)
local.plot.field(corr, mesh4, zlim=c(0.1, 1))
points(6.90, 10.21)
plot(poly.water, add=T)


## --------------------------------------------------------------------------------
## Remove this and update the code to use the new functionality in INLA
# - Use inla.barrier.pcmaterns instead of the code below
source("functions-barriers-dt-models-march2017.R")


## --------------------------------------------------------------------------------
mesh = mesh4
tl = length(mesh$graph$tv[,1])
# - the number of triangles in the mesh
posTri = matrix(0, tl, 2)
for (t in 1:tl){
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}
posTri = SpatialPoints(posTri)
# - The positions of the triangles
water = over(poly.water, SpatialPoints(posTri), returnList=T)
# - checking which mesh triangles are inside the normal area
water = unlist(water)
Omega = list(water, setdiff(1:tl, water))
Omega.SP = dt.polygon.omega(mesh, Omega)
# - creates polygons for the different areas
# - - the first is water/normal area 
# - - the second is Barrier area
str(Omega.SP, 1)


## --------------------------------------------------------------------------------
corr = local.find.correlation(Q, loc = c(7,10.3), mesh)
local.plot.field(corr, mesh, zlim=c(0.1, 1))
points(6.90, 10.21)
plot(Omega.SP[[2]], add=T, col="grey")
# - adding the land polygon filled with grey
# - this hides the field on land


## --------------------------------------------------------------------------------
local.plot.field(corr, mesh, xlim = c(5, 9), ylim = c(8, 12), zlim=c(0.1, 1))
points(6.90, 10.21)
plot(Omega.SP[[2]], add=T, col="grey")


## --------------------------------------------------------------------------------
Q.function = dt.create.Q(mesh, Omega, fixed.ranges = c(NA, 0.5))
# - the 0.5-fixed range is for the barrier area
# - - it is not sensitive to the exact value here, 
#     just make it "small"


## --------------------------------------------------------------------------------
r = 3
# - some chosen range (in the water area)
sigma = 1
# - some chosen sigma scaling

Q = Q.function(theta = c(log(sigma), log(r)))
# - the precision matrix for fixed ranges
# - Q is a function of the hyperparameters theta = c( log(sigma), log(range1), log(range2),...)


## ---- results="hold"-------------------------------------------------------------
Q = Q.function(theta = c(log(1), log(r)))
sd = diag(inla.qinv(Q))
local.plot.field(sd, mesh)
plot(Omega.SP[[2]], add=T, col="grey")
# - we only care about our study area


## --------------------------------------------------------------------------------
r = 4
Q = Q.function(theta = c(log(1), log(r)))
corr = local.find.correlation(Q, loc = c(5,5), mesh)
local.plot.field(corr, mesh, zlim=c(0.13, 1))
plot(Omega.SP[[2]], add=T, col="grey")


## --------------------------------------------------------------------------------
r = 3
Q = Q.function(theta = c(log(1), log(r)))
corr = local.find.correlation(Q, loc = c(5,9), mesh)
local.plot.field(corr, mesh, zlim=c(0.13, 1))
plot(Omega.SP[[2]], add=T, col="grey")

