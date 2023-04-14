## ----setup, include=FALSE-----------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- message= F--------------------------------------------------
library(INLA)
library(fields)
library(rgeos)
library(ggplot2)


## -----------------------------------------------------------------
## Load data
## After you have downloaded/saved the file
load(file = "data/WebSiteData-Archipelago.RData")

## What is loaded
# - poly.water is our study area
# - df is our dataframe to be analysed
# - dat is the orginial dataframe
str(poly.water, 1)


## -----------------------------------------------------------------
df$y = df$y.smelt


## -----------------------------------------------------------------
## Set the max length of triangles edge
max.edge = 0.6
## Set the length of the boundary extension
bound.outer = 4.6
## Build the mesh
mesh = inla.mesh.2d(boundary = poly.water,
                    loc=cbind(df$locx, df$locy),
                    max.edge = c(1,5)*max.edge,
                    cutoff = 0.06,
                    offset = c(max.edge, bound.outer))


## -----------------------------------------------------------------
plot(mesh, lwd=0.5) 
points(df$locx, df$locy, col="red")


## -----------------------------------------------------------------
## Locations for observations
locations = cbind(df$locx, df$locy)
## Projection matrix
A.i.s = inla.spde.make.A(mesh, loc = locations)


## ---- results='hold'----------------------------------------------
cat('Dimension of A: ', dim(A.i.s), '\n')

cat('Number of mesh points: ', mesh$n, '\n')

cat('Number of locations: ', dim(locations)[1], '\n')


## -----------------------------------------------------------------
stk <- inla.stack(data=list(y=df$y, e=df$exposure), # data and offset
                  effects=list(s= 1:mesh$n, # spatial random effect
                               iidx=1:nrow(df), # iid random effect
                               data.frame(m=1, df[ ,5:11])), #covariates
                  A=list(A.i.s, 1, 1), # projection matrix
                  remove.unused = FALSE, tag='est')



## -----------------------------------------------------------------
## Overall list
M = list()

# Each model is a list itself
M[[1]] = list()
M[[1]]$shortname = "stationary-no-cov"
M[[2]] = list()
M[[2]]$shortname = "stationary-all-cov"
M[[3]] = list()
M[[3]]$shortname = "barrier-no-cov"
M[[4]] = list()
M[[4]]$shortname = "barrier-all-cov"


## -----------------------------------------------------------------
## Prior for range uses half the study area, which is
## ... approximately 0.5*diff(range(df$locy)) = 5.8435
spde = inla.spde2.pcmatern(mesh, prior.range = c(6, .5), prior.sigma = c(3, 0.01))


## -----------------------------------------------------------------
hyper.iid = list(prec = list(prior = 'pc.prec', param = c(3, 0.01))) 


## -----------------------------------------------------------------
## No covariates
M[[1]]$formula = y ~ -1 + m + 
  f(s, model=spde) +                      # spatial random effect
  f(iidx, model="iid", hyper=hyper.iid)   # iid random effect


## -----------------------------------------------------------------
## Covariates
M[[2]]$formula = as.formula(paste( "y ~ -1 + ",paste(colnames(df)[5:11], collapse = " + ")))

## Add random effect
M[[2]]$formula = update(M[[2]]$formula, .~. + m + 
                          f(s, model=spde) + 
                          f(iidx, model="iid", hyper=hyper.iid))



## -----------------------------------------------------------------
print(M[[1]])


## -----------------------------------------------------------------
print(M[[2]])


## -----------------------------------------------------------------
# Number of triangles of the mesh
tl = length(mesh$graph$tv[,1])

## Initialize a matrix containing the central coordinates of each triangle's 
posTri = matrix(0, tl, 2)

for(t in 1:tl){
  # Take the vertex of triangles
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  # Compute center of each triangle
  posTri[t,] = colMeans(temp)[c(1,2)]
}

## Convert to Spatial Points
posTri = SpatialPoints(posTri)

## Intersection between mesh points and sea points contained in 
## ... poly.water
normal = over(poly.water, posTri, returnList=T)

## Remove the sea triangles from all triangles to obtain the barrier ones
normal = unlist(normal)
barrier.triangles = setdiff(1:tl, normal)

## Build a barrier, this obj contains all the polygons composing the islands
## This will be used just for plotting
poly.barrier = inla.barrier.polygon(mesh, barrier.triangles)



## -----------------------------------------------------------------
barrier.model = inla.barrier.pcmatern(mesh, barrier.triangles = 
                                        barrier.triangles, 
                                      prior.range = c(6, .5), 
                                      prior.sigma = c(3, 0.01))


## -----------------------------------------------------------------
## No Covariates
M[[3]]$formula = y~ -1 + m + 
  f(s, model=barrier.model) + 
  f(iidx, model="iid", hyper=hyper.iid)

## With Covariates
## Only covariates
M[[4]]$formula = as.formula(paste( "y ~ -1 + ",paste(colnames(df)[5:11], collapse = " + ")))
## Add random effects
M[[4]]$formula = update(M[[4]]$formula, .~. +m + 
                          f(s, model=barrier.model) + 
                          f(iidx, model="iid", hyper=hyper.iid))


## -----------------------------------------------------------------
print(M[[3]])
print(M[[4]])


## -----------------------------------------------------------------
local.find.correlation = function(Q, location, mesh) {
  ## Vector of standard deviations
  sd = sqrt(diag(inla.qinv(Q)))
  
  ## Create a fake A matrix, to extract the closest mesh node index
  A.tmp = inla.spde.make.A(mesh=mesh, 
                           loc = matrix(c(location[1],location[2]),1,2))
  
  ## Index of the closest node
  id.node = which.max(A.tmp[1, ])
  
  
  print(paste('The location used was c(', 
              round(mesh$loc[id.node, 1], 4), ', ', 
              round(mesh$loc[id.node, 2], 4), ')' ))
  
  ## Solve a matrix system to find the column of the covariance matrix
  Inode = rep(0, dim(Q)[1]) 
  Inode[id.node] = 1
  covar.column = solve(Q, Inode)
  # compute correaltions
  corr = drop(matrix(covar.column)) / (sd*sd[id.node])
  return(corr)
}



## -----------------------------------------------------------------
local.plot.field = function(field, mesh, xlim, ylim, ...){
  # Error when using the wrong mesh
  stopifnot(length(field) == mesh$n)
  
  # Choose plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim = poly.water@bbox[1, ] 
  if (missing(ylim)) ylim = poly.water@bbox[2, ]
  
  # Project the mesh onto a 300x300 grid
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  
  # Do the projection 
  field.proj = inla.mesh.project(proj, field)
  
  # Plot it
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, ...)  
}



## -----------------------------------------------------------------
# theta = c(log(range), log(sigma))
Q = inla.spde2.precision(spde, theta = c(log(6),log(3)))

corr = local.find.correlation(Q, loc = c(8,10), mesh)

local.plot.field(corr, mesh, zlim=c(0.1, 1))

points(8.01, 9.97)
plot(poly.barrier, add=T, col='grey', main = 'Stationary Model')
title(main = 'Stationary Model')


## -----------------------------------------------------------------
Q = inla.rgeneric.q(barrier.model, "Q", theta = c(0, log(6)))

corr = local.find.correlation(Q, loc = c(8,10), mesh)

local.plot.field(corr, mesh, zlim=c(0.1, 1))
points(8.01, 9.97)
plot(poly.barrier, add=T, col='grey', main = 'Barrier Model')
title(main = 'Barrier Model')


## -----------------------------------------------------------------
M[[1]]$init = c(2.509,1.199,-0.574)
M[[2]]$init = c(1.162,0.313,-0.627)
M[[3]]$init = c(0.833,2.244,-0.471)
M[[4]]$init = c(0.044,1.274,-0.596)


## ---- warning = F-------------------------------------------------
for (i in 1:length(M)){
  print(paste("Running:  ", M[[i]]$shortname))
  M[[i]]$res = inla(M[[i]]$formula,
                    data=inla.stack.data(stk),
                    control.predictor=list(A=inla.stack.A(stk)),
                    family="poisson", E = e,
                    control.inla= list(int.strategy = "eb"),
                    control.mode=list(restart=T, theta=M[[i]]$init))  
}


## -----------------------------------------------------------------
for (i in 1:length(M)){
  print(paste(round(M[[i]]$res$internal.summary.hyperpar$mode, 3), collapse = ','))
}


## -----------------------------------------------------------------
# Set up dataframe with relevant results
res = M[[2]]$res$summary.fixed[ ,c(4,3,5)]
res = rbind(res, M[[4]]$res$summary.fixed[ ,c(4,3,5)])

# Change the name of the column:
# E = Estimates
# L = Lower
# U = Upper
colnames(res) = c("E", "L", "U")
rownames(res)=NULL

# number of covariates
n.covar = nrow(M[[2]]$res$summary.fixed)

# specify a factor corresponding to the model:
# MS = stationary
# MB = barrier
res$model = factor(rep(c("MS", "MB"), each=n.covar, 
                       levels = c("MS", "MB")))

# specify a factor corresponding to the covariates
res$covar = factor(rep(rownames(M[[2]]$res$summary.fixed), 2))

# Plot it
ggplot(res, aes(x = model, y = E)) +
  facet_wrap(~covar, scales = "free_y") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  xlab(NULL) + ylab(NULL) 


## -----------------------------------------------------------------
## Quantiles considered
quantile.words = c(".5", ".025", ".975", ".1", ".9")


local.make.quantile.table = function(stat, barr, quantile.words){
  
  ## Initialize the dataframe
  names = names(stat_nocov)
  pos.hyp = data.frame(name = rep(names,2))
  qnames = paste0("q", quantile.words)
  for (name in qnames) {
    pos.hyp[[name]] = NA
  }
  
  ## Fill the quantiles information
  for(i in 1:3){
    pos.hyp[i,2:6] = inla.qmarginal(as.numeric(quantile.words), stat[[i]])
    pos.hyp[(i+3),2:6] = inla.qmarginal(as.numeric(quantile.words), barr[[i]])
  }
  
  ## Add a factor representing the model
  pos.hyp$model = factor(rep(c("MS", "MB"), each=3, 
                             levels = c("MS", "MB")))
  return(pos.hyp)
}



## -----------------------------------------------------------------
## Results stationary model
stat_nocov = list(
  sigma.epsilon = inla.tmarginal(function(x) exp(-0.5*x),
                                 M[[1]]$res$internal.marginals.hyperpar[[3]]),
  sigma.u = inla.tmarginal(function(x) x,
                           M[[1]]$res$marginals.hyperpar[[2]]),
  range = inla.tmarginal(function(x) x,
                         M[[1]]$res$marginals.hyperpar[[1]]))

## Results barrier model
barr_nocov = post = list(
  sigma.epsilon = inla.tmarginal(function(x) exp(-0.5*x), 
                                 M[[3]]$res$internal.marginals.hyperpar[[3]]),
  sigma.u = inla.tmarginal(function(x) exp(x), 
                           M[[3]]$res$internal.marginals.hyperpar[[1]]),
  range = inla.tmarginal(function(x) exp(x), 
                         M[[3]]$res$internal.marginals.hyperpar[[2]]))



## ---- fig.height=4------------------------------------------------
pos.hyp_nc = local.make.quantile.table(stat_nocov, barr_nocov, quantile.words)

ggplot(pos.hyp_nc, aes(x = model, y = q.5)) +
  facet_wrap(~name, scales = "free_y") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymax = q.975, ymin = q.025)) +
  geom_errorbar(aes(ymax = q.9, ymin = q.1), col="black") +
  xlab(NULL) + ylab(NULL) 




## ---- fig.height=4------------------------------------------------
## Results stationary model
stat_cov = list(
  sigma.epsilon = inla.tmarginal(function(x) exp(-0.5*x),
                                 M[[2]]$res$internal.marginals.hyperpar[[3]]),
  sigma.u = inla.tmarginal(function(x) x,
                           M[[2]]$res$marginals.hyperpar[[2]]),
  range = inla.tmarginal(function(x) x,
                         M[[2]]$res$marginals.hyperpar[[1]]))

## Results barrier model
barr_cov = post = list(
  sigma.epsilon = inla.tmarginal(function(x) exp(-0.5*x), 
                                 M[[4]]$res$internal.marginals.hyperpar[[3]]),
  sigma.u = inla.tmarginal(function(x) exp(x), 
                           M[[4]]$res$internal.marginals.hyperpar[[1]]),
  range = inla.tmarginal(function(x) exp(x), 
                         M[[4]]$res$internal.marginals.hyperpar[[2]]))


pos.hyp_c = local.make.quantile.table(stat_cov, barr_cov, quantile.words)

ggplot(pos.hyp_c, aes(x = model, y = q.5)) +
  facet_wrap(~name, scales = "free_y") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymax = q.975, ymin = q.025)) +
  geom_errorbar(aes(ymax = q.9, ymin = q.1), col="black") +
  xlab(NULL) + ylab(NULL) 




