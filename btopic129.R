## ----setup, include=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- message = F, warning=FALSE-------------------------------
library(INLA)
library(fields)
library(rgeos)
library(sp)
library(colorRamps)


## --------------------------------------------------------------
## Load data
load(file = "data/WebSiteData-Archipelago.RData")
# - if you have saved the file locally

## What is loaded
# - poly.water is our study area
# - df is our dataframe to be analysed
# - dat is the orginial dataframe
str(poly.water, 1)


## --------------------------------------------------------------
# Find the boundaries of the reduced study area
xmax = poly.water@bbox[1,2]
xmin = xmax/2
ymax = poly.water@bbox[2,2]

# Define the reduced study area as a Polygon 
study_area = Polygon(cbind(c(0, xmin, xmin, 0, 0),
                           c(0, 0, ymax, ymax, 0)), hole = F)
study_area = SpatialPolygons(list(Polygons(list(study_area),'1')))

# Intersect with the entire study area
shape <- gBuffer(poly.water, byid=TRUE, width=0) #for gIntersection operation
shape2 <- gIntersection(shape, study_area)
poly.water = shape2


# Take a look
plot(poly.water, main = 'Reduced study area')


## ---- fig.width=8, fig.height=8--------------------------------
set.seed(2018)
set.inla.seed = 2018

max.edge = 0.6
bound.outer = 3.6
prmesh1 = inla.mesh.2d(boundary = poly.water,
 #                      loc=cbind(df$locx, df$locy),
                       max.edge = c(1,5)*max.edge,
                       cutoff = 0.06,
                       offset = c(max.edge, bound.outer))

# take a look
plot(prmesh1, lwd = 0.7, main="Triangulation of the study area")
title(main = "Triangulation of the study area")
#points(df$locx, df$locy, pch = 21, bg = 2)



## --------------------------------------------------------------
# Total number of triangles 
tl = length(prmesh1$graph$tv[,1])
# Index vector of all triangles
all.triangles = 1:tl

# Matrix containing the location of the barycenters
barycenters = matrix(0, nrow = tl, ncol = 2)
for(t in 1:tl){
  temp = prmesh1$loc[prmesh1$graph$tv[t, ], ]
  barycenters[t,] = colMeans(temp)[c(1,2)]
}

# Transform it in SpatialPoints to intersect it with the barrier easily
barycenters = SpatialPoints(barycenters)

# Intersect
not.barrier = over(poly.water, barycenters, returnList=T)
not.barrier = unlist(not.barrier)

# Difference
barrier.triangles = setdiff(all.triangles, not.barrier)

# This object will be usefull later, just to plot.
poly.barrier = inla.barrier.polygon(prmesh1, barrier.triangles)


## --------------------------------------------------------------
loc.data <- spsample(poly.water, n = 1000, type = "random")
coords <- loc.data@coords


## --------------------------------------------------------------
k <- 4 
range <- sqrt(8)
sigma <- 1
rho <- 0.7

# input parmaters for the INLA function
theta <- c(log(range),log(sigma))
# mesh for the simulation
simulation.mesh <- prmesh1


## ---- warning=FALSE--------------------------------------------
# Create spde to simulate
# Ignore range and sigma priors, the function doesn't work without, you can put
# there everything
simulation.spde <- inla.barrier.pcmatern(mesh=simulation.mesh,
                                         barrier.triangles = barrier.triangles,
                                         prior.range = c(0.1,0.1), 
                                         prior.sigma = c(0.1,0.1))
  
# Create precision  
simulation.Q <- inla.rgeneric.q(simulation.spde, "Q", theta=theta)

# Create projection matrix
simulation.A <- inla.mesh.project(mesh=simulation.mesh, loc=coords)$A


## --------------------------------------------------------------
# simulate 4 independent samples
x <- inla.qsample(k, simulation.Q, seed=0, constr=simulation.spde$f$extraconstr)

# initialize a matrix
x.hat = matrix(NA, nrow = nrow(simulation.A), ncol = ncol(x))

# project each sample on the data locations
for (j in 1:ncol(x)){
  x.hat[, j] <- drop(simulation.A%*%x[,j])
}

# aggregate to create an AR(1)
for (j in 2:k){
  x.hat[,j] <- rho*x.hat[,j-1] + sqrt(1-rho^2)*x.hat[,j]
}

# Take a look

rbPal <- colorRampPalette(matlab.like2(100))
par(mfrow=c(2,2), mar=c(0,1,1,0))
for (j in 1:k){
  Col <- rbPal(100)[as.numeric(cut(x.hat[,j], breaks = 100))]
  plot(coords, pch = 20, col = Col, xlab = '', ylab = '', 
       main = paste0('t = ', j), axes = F)
  plot(poly.barrier, add = T, col = 'grey')
}



## --------------------------------------------------------------
set.seed(2)

n <- nrow(coords)
# sample the categorical covariate
w <- sample(-1:1, n*k, replace=TRUE)
# take a look
table(w)


sd.y <- 0.1
# generate the target variable
y <- x.hat + w + rnorm(n*k, 0, sd.y) 

w <- factor(w)
#see the mean  of y which have the same value of w
tapply(y, w, mean)



## --------------------------------------------------------------
isel <- sample(1:(n*k), n*k/2)
dat <- data.frame(y=as.vector(y),
                  w=w,
                  time=rep(1:k, each=n), #notice the time
                  xcoo=rep(coords[,1], k),
                  ycoo=rep(coords[,2], k))[isel, ]



## --------------------------------------------------------------
spde <- inla.barrier.pcmatern(mesh=prmesh1,
                              barrier.triangles = barrier.triangles,
                              prior.range=c(0.5, 0.01), #P(range < 0.5) = 0.01
                              prior.sigma=c(1, 0.01)) #P(sigma > 1) = 0.01

# fix the other two priors
rho.prior <- list(theta=list(prior='pccor1', param=c(0, 0.9))) #P(rho > 0) = 0.9
theta.prior <- list(prior='pc.prec', param=c(1, 0.01)) # P(tau < 1) = 0.01

# create an index to pass to the stack and to f() later
iset <- inla.spde.make.index('spatio.temp', n.spde = prmesh1$n, n.group = k)



## --------------------------------------------------------------
# Create projection matrix (it includes also the time!)
A <- inla.spde.make.A(mesh=prmesh1, 
                      loc=cbind(dat$xcoo, dat$ycoo),
                      group=dat$time)

# Build the stack
stack <- inla.stack(tag = 'stdata', 
                   data = list(y=dat$y), # target variable 
                   A = list(A,1), # projection matrices, 1 = no projection
                   effects = list(iset, w = dat$w)) # list of effects



## --------------------------------------------------------------
formula <- y ~ 0 + w +f(spatio.temp, model=spde, group=spatio.temp.group,
                        control.group=list(model='ar1', hyper=rho.prior))


## ---- warning=F, message=FALSE---------------------------------
init = c(4.689,  0.977, -0.007,  1.714)
res <- inla(formula, data=inla.stack.data(stack),
            control.predictor=list(compute=TRUE, A=inla.stack.A(stack)), 
            control.family=list(hyper=list(theta=theta.prior)), 
            control.fixed=list(expand.factor.strategy='inla'),
            control.compute = list(return.marginals=F, config=T),
            control.inla= list(int.strategy = "eb"), # line to add to speed up
            control.mode=list(restart=T, theta=init)) # line to add to speed up

res$internal.summary.hyperpar$mode


## --------------------------------------------------------------
tapply(dat$y, dat$w, mean)  # Observed mean for each covariate level
round(res$summary.fixed,4) # Coefficients of the categorical variable


## --------------------------------------------------------------
names(res$marginals.hyper)


## --------------------------------------------------------------
par(mfrow=c(2,2), mar=c(3,3,1,0.1), mgp=2:0)
plot(res$marginals.hyper[[1]], type='l', 
     xlab=names(res$marginals.hyper)[1], ylab='Density',xlim=c(0,200))
abline(v=c(1/sd.y^2, log(range),
           log(sigma), rho)[1], col=2)
for (j in 2:4) {
  plot(res$marginals.hyper[[j]], type='l', 
       xlab=names(res$marginals.hyper)[j], ylab='Density')
  abline(v=c(1/sd.y^2, log(range),
             log(sigma), rho)[j], col=2) 
}


## --------------------------------------------------------------
names(res$summary.random$spatio.temp)


## --------------------------------------------------------------
projgrid <- inla.mesh.projector(prmesh1, 
                                xlim=range(coords[,1]),
                                ylim=range(coords[,2]), 
                                dims=c(100,100)) #produce the projection grid

#produce the prediction result for each time
xmean <- list()
for (j in 1:k){
  xmean[[j]] <- inla.mesh.project(projgrid,
                        res$summary.random$spatio.temp$mean[iset$spatio.temp.group==j])
} 
  
#plot the prediction result

par(mfrow=c(2,2), mar = c(1,1,1,1))
for (i in 1:k){
  image.plot(list(x = projgrid$x, y=projgrid$y, z = xmean[[i]]), xlim=c(0,10),
             ylim=c(0,13), main = paste0('time = ', i), axes = F, zlim=c(-10, 10))
  plot(poly.barrier, add=T, col='grey')
}



## --------------------------------------------------------------
par(mfrow=c(1,2), mar=c(1,1,1,1))

# plot true one
rbPal <- colorRampPalette(matlab.like2(100))
Col <- rbPal(100)[as.numeric(cut(x.hat[,1], breaks = 100))]
plot(coords, pch = 20, col = Col, xlab = '', ylab = '', cex = 1, 
     axes = F, main = 'True')
plot(poly.barrier, add = T, col  ='grey')

# plot estimate
image.plot(list(x = projgrid$x, y=projgrid$y, z = xmean[[1]]), 
           xlim=c(0,10), ylim=c(0,13), 
           col = matlab.like2(100), axes = F, main = 'Posterior mean')
plot(poly.barrier, add=T, col='grey', xlim = c(0,10))


