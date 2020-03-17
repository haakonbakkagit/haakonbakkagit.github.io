## ----setup, include=FALSE-------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE----------------------------------
library(INLA)
library(fields)
library(rgdal)
library(viridisLite)

set.seed(2016)
set.inla.seed = 2016



## -------------------------------------------------------------------
smalldist = 0.2
# - the width of the opening in the barrier
width = 0.5
# - The width/thickness of the barrier

local.square.polygon = function(xlim, ylim){
# - output is a square
  xlim = range(xlim); ylim = range(ylim)
  corner1 = c(xlim[1], ylim[2])
  corner2 = c(xlim[2], ylim[1])
  poly = Polygon(rbind(corner1, c(corner1[1], corner2[2]), corner2, c(corner2[1], corner1[2]), corner1), hole = FALSE)
  return(SpatialPolygons(list(Polygons(list(poly), ID = runif(1)))))
}

poly1 = local.square.polygon(xlim=c(-1, 5-smalldist/2), 
                          ylim=5+width*c(-.5, .5))
poly2 = local.square.polygon(xlim=c(5+smalldist/2, 11), 
                          ylim=5+width*c(-.5, .5))
poly.original = SpatialPolygons(c(poly1@polygons, poly2@polygons))

plot(poly.original, main="Barrier area polygon")


## -------------------------------------------------------------------
max.edge.length = 0.4
# - The coarseness of the finite element approximation
# - Corresponds to grid-square width in discretisations
# - - Except that finite element approximations are better
# - Should be compared to size of study area
# - Should be less than a fourth of the estimated (posterior) 
#   spatial range
# - Up to 8x computational time when you halve this value



## -------------------------------------------------------------------
loc1 = matrix(c(0,0, 10,0, 0,10, 10,10), 4, 2, byrow = T)
# - This defines the extent of the interior part of the mesh
# - In an application, if you want the mesh to depend on your 
#   data locations, you may use those locations instead
seg = inla.sp2segment(poly.original)
# - Transforms a SpatialPolygon to an "inla polygon"
mesh = inla.mesh.2d(loc=loc1, interior = seg, 
                    max.e = max.edge.length, offset=1)
# - The INLA mesh constructor, used for any INLA-SPDE model



## -------------------------------------------------------------------
tl = length(mesh$graph$tv[,1])
# - the number of triangles in the mesh
posTri = matrix(0, tl, 2)
for (t in 1:tl){
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}
posTri = SpatialPoints(posTri)
# - the positions of the triangle centres

barrier = over(poly.original, SpatialPoints(posTri), returnList=T)
# - checking which mesh triangles are inside the barrier area
barrier = unlist(barrier)
poly.barrier = inla.barrier.polygon(mesh, barrier.triangles = barrier)
# - the Barrier model's polygon
# - in most cases this should be the same as poly.original


## -------------------------------------------------------------------
plot(mesh, main="Mesh and Omega")
plot(poly.barrier, add=T, col='lightblue')
plot(mesh, add=T)
points(loc1)


## -------------------------------------------------------------------
local.plot.field = function(field, ...){
  xlim = c(2, 8); ylim = xlim;
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 grid 
  #   for plots
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, col = plasma(17), ...)  
  # - Use image.plot to get nice colors and legend
}
print(mesh$n)
# - This is the appropriate length of the field variable


## -------------------------------------------------------------------
barrier.model = inla.barrier.pcmatern(mesh, barrier.triangles = barrier, prior.range = c(1.44, 0.5), prior.sigma = c(0.7, 0.5), range.fraction = 0.1)
# - Set up the inla model, including the matrices for solving the SPDE


## -------------------------------------------------------------------
range = 3
# - the spatial range parameter
Q = inla.rgeneric.q(barrier.model, "Q", theta = c(0, log(range)))
# - the precision matrix for fixed ranges
# - Q is a function of the hyperparameters theta = c( log(sigma), log(range1), log(range2),...)


## -------------------------------------------------------------------
u = inla.qsample(n=1, Q=Q, seed = set.inla.seed)
u = u[ ,1]
# - access the first sample

local.plot.field(u, main="The true (simulated) spatial field")
plot(poly.barrier, add=T, col='grey')
# - Overlay the barrier


## -------------------------------------------------------------------
num.try = 500 
# - try to sample this number of data locations
loc.try = matrix(runif(num.try*2, min=2, max=8), 
                         num.try, 2)
# - locations sampled inside the barrier will be removed 
#   in a few lines
temp = SpatialPoints(loc.try)
loc.ok = is.na(over(temp, poly.barrier))
# - only allow locations that are not inside the Barrier area
loc.data = loc.try[loc.ok, ]
A.data = inla.spde.make.A(mesh, loc.data)
# - the projector matrix required for any spatial model
# - this matrix can transform the field-defined-on-the-mesh 
#   to the field-defined-on-the-data-locations
c(dim(A.data), mesh$n, nrow(loc.data))
# - shows that the dimensions are correct
u.data = A.data %*% u
# - project the field from the finite element  
#   representation to the data locations


## -------------------------------------------------------------------
df = data.frame(loc.data)
# - df is the dataframe used for modeling
names(df) = c('locx', 'locy')
sigma.u = 1
# - size of the random effect
# - feel free to change this value
sigma.epsilon = 0.2
# - size of the iid noise in the Gaussian likelihood
# - feel free to change this value
df$y = drop(sigma.u*u.data + sigma.epsilon*rnorm(nrow(df)))
# - sample observations with gaussian noise


## -------------------------------------------------------------------
summary(df)


## -------------------------------------------------------------------
stk <- inla.stack(data=list(y=df$y), A=list(A.data, 1),
                  effects=list(s=1:mesh$n, 
                               intercept=rep(1, nrow(df))), 
                  remove.unused = FALSE, tag='est')
# - this is the common stack used in INLA SPDE models
# - see the SPDE-tutorial
# - - http://www.r-inla.org/examples/tutorials/spde-tutorial


## -------------------------------------------------------------------
model.stat = inla.spde2.pcmatern(mesh, prior.range = c(1, 0.5), prior.sigma = c(1, 0.5))
# - Set up the model component for the spatial SPDE model: 
#   Stationary Matern model
# - I assume you are somewhat familiar with this model

formula <- y ~ 0 + intercept + f(s, model=model.stat)
# - Remove the default intercept
# - - Having it in the stack instead improves the numerical 
#     accuracy of the INLA algorithm
# - Fixed effects + random effects

res.stationary <- inla(formula, data=inla.stack.data(stk),
            control.predictor=list(A = inla.stack.A(stk)),
            family = 'gaussian',
            control.family = list(hyper = list(prec = list(
              prior = "pc.prec", fixed = FALSE, 
              param = c(0.2,0.5)))),
            control.mode=list(restart=T, theta=c(4,-1.7,0.25)))


## -------------------------------------------------------------------
summary(res.stationary)


## -------------------------------------------------------------------
local.plot.field(res.stationary$summary.random$s$mean,
          main="Spatial estimate with the stationary model")
# - plot the posterior spatial marginal means
# - we call this the spatial estimate, or the smoothed data
plot(poly.barrier, add=T, col='grey')
# - Posterior spatial estimate using the stationary model


## -------------------------------------------------------------------
formula2 <- y ~ 0 + intercept + f(s, model=barrier.model)
# - The spatial model component is different from before
# - The rest of the model setup is the same as in the stationary case!
# - - e.g. the inla(...) call below is the same, 
#     only this formula is different


## -------------------------------------------------------------------
res.barrier = inla(formula2, data=inla.stack.data(stk),
       control.predictor=list(A = inla.stack.A(stk)),
       family = 'gaussian',
       control.family = list(hyper = list(prec = list(
             prior = "pc.prec", fixed = FALSE, 
             param = c(0.2,0.5)))),
       control.mode=list(restart=T, theta=c(3.2, 0.4, 1.6)))


## -------------------------------------------------------------------
summary(res.barrier)


## -------------------------------------------------------------------
local.plot.field(res.barrier$summary.random$s$mean, 
                 main="Spatial posterior for Barrier model")
# - plot the posterior spatial marginal means
# - we call this the spatial (smoothing) estimate
plot(poly.barrier, add=T, col='grey')
# - Posterior spatial estimate using the Barrier model


## -------------------------------------------------------------------
res.barrier$summary.hyperpar


## -------------------------------------------------------------------
tmp = inla.tmarginal(function(x) exp(x), res.barrier$marginals.hyperpar[[2]]) 
plot(tmp, type = "l", xlab = "sigma", ylab = "Density")
xvals = seq(0, 10, length.out=1000)
lambda = 0.99; lines(xvals, 3*exp(-lambda*xvals), lty='dashed')
abline(v=1, col="blue")


## -------------------------------------------------------------------
tmp = inla.tmarginal(function(x) exp(x), res.barrier$marginals.hyperpar[[3]]) 
plot(tmp, type = "l", xlab = "r", ylab = "Density")
xvals = seq(0, 10, length.out=1000)
lambda = 1.00; lines(xvals, 3*exp(-lambda*xvals), lty='dashed')
abline(v=range, col="blue")

