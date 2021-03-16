## ----setup, include=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## --------------------------------------------------------------
set.seed(1)
n = 20
x_loc = cbind(runif(n), runif(n))


## ---- include=F------------------------------------------------
library(INLA)


## --------------------------------------------------------------
mesh1 = inla.mesh.2d(loc = x_loc, max.edge = 0.1)
mesh2 = inla.mesh.2d(loc = x_loc, max.edge = c(0.1, 0.4))
mesh3 = inla.mesh.2d(loc.domain = x_loc, max.edge = c(0.1, 0.4))

## ---- echo = F, fig.height=3-----------------------------------
par(mfrow = c(1,3))
plot(mesh1, main = 'mesh1')
title(main = 'mesh1')
points(x_loc, pch = 21, bg = 2)
plot(mesh2, main = 'mesh2')
title(main = 'mesh2')
points(x_loc, pch = 21, bg = 2)
plot(mesh3, main = 'mesh3')
title(main = 'mesh3')
points(x_loc, pch = 21, bg = 2)



## --------------------------------------------------------------
x_loc_reduce = x_loc[1:5,]
m1 <- inla.mesh.2d(x_loc_reduce, max.edge=c(0.5, 0.5))
m2 <- inla.mesh.2d(x_loc_reduce, max.edge=c(0.5, 0.5), cutoff = 0.1)

m3 <- inla.mesh.2d(x_loc_reduce, max.edge=c(0.1, 0.5), cutoff = 0.1)
m4 <- inla.mesh.2d(x_loc_reduce, max.edge=c(0.1, 0.5), cutoff = 0.1,
                   offset=c(0, 0.3))
m5 <- inla.mesh.2d(x_loc_reduce, max.edge=c(0.1, 0.5), cutoff = 0.1,
                   offset=c(0.2, 0.3))
m6 <- inla.mesh.2d(x_loc_reduce, max.edge=c(0.1, 0.5), cutoff = 0.1,
                   offset=c(0.2, 0.7), min.angle = 10)



## ---- echo = F, fig.height=5-----------------------------------
par(mfrow = c(2,3)) 
for (i in 1:6) {
  plot(get(paste('m', i, sep='')), main = paste('m',i,sep=''))
  title(main = paste('m',i,sep=''))
  points(x_loc_reduce, pch = 21, bg = 2)
}


## --------------------------------------------------------------
boundary = cbind(c(0.2,0.2,1,1), c(0.1,1,1,0.1))

m7 <- inla.mesh.2d(, boundary, max.edge=c(0.3, 0.5), 
                   offset=c(0.03, 0.5), cutoff=0.1)
m8 <- inla.mesh.2d(, boundary, max.edge=c(0.3, 0.5), n=5, 
                   offset=c(.05,.1))
m9 <- inla.mesh.2d(, boundary, max.edge=c(.3, 0.5), n=7, 
                   offset=c(.01,.3))
m10 <- inla.mesh.2d(, boundary, max.edge=c(.3, 0.5), n=4, 
                   offset=c(.05,.3))



## ---- echo = F, fig.height=6-----------------------------------
boundary = cbind(c(0.2,0.2,1,1,0.2), c(0.1,1,1,0.1,0.1))

par(mfrow = c(2,2)) 
for (i in 7:10) {
  plot(get(paste('m', i, sep='')), main = paste('m',i,sep=''))
  title(main = paste('m',i,sep=''))
  points(boundary, type = 'l', col = 3, lwd = 2)
  points(x_loc_reduce, pch = 21, bg = 2)
}



## --------------------------------------------------------------
bnd9 = inla.nonconvex.hull(x_loc, convex = 0.05)
bnd10 = inla.nonconvex.hull(x_loc, convex = 0.09)
bnd11 = inla.nonconvex.hull(x_loc, convex = 0.2)

m9 = inla.mesh.2d(boundary = bnd9, max.edge = 0.05)
m10 = inla.mesh.2d(boundary = bnd10, max.edge = 0.05)
m11 = inla.mesh.2d(boundary = bnd11, max.edge = 0.05)


## ---- echo = F, fig.height=3-----------------------------------
par(mfrow = c(1,3))
lines(m9$segm$bnd, m9$loc, add = F)
points(x_loc, pch = 21, bg = 2)
lines(m10$segm$bnd, m10$loc, add = F)
points(x_loc, pch = 21, bg = 2)
lines(m11$segm$bnd, m11$loc, add = F)
points(x_loc, pch = 21, bg = 2)


## --------------------------------------------------------------
m9 = inla.mesh.2d(boundary = bnd9, max.edge = c(0.09, 0.4), cutoff = c(0.01, 0.01))
m10 = inla.mesh.2d(boundary = bnd10, max.edge = c(0.09, 0.4), 
                   offset = c(0, -0.2))
m11 = inla.mesh.2d(boundary = bnd11, max.edge = c(0.09, 0.4), 
                   min.angle = 0.05)


## ---- echo = F, fig.height=3-----------------------------------
par(mfrow = c(1,3))
plot(m9, main = 'm9'); title(main = 'm9')
points(x_loc, pch = 21, bg = 2)
plot(m10, main = 'm9'); title(main = 'm10')
points(x_loc, pch = 21, bg = 2)
plot(m11, main = 'm9'); title(main = 'm11')
points(x_loc, pch = 21, bg = 2)


## --------------------------------------------------------------
globe1 = inla.mesh.create(globe = 1)
globe2 = inla.mesh.create(globe = 4)
globe3 = inla.mesh.create(globe = 10)


## ---- echo = F, fig.height=4-----------------------------------
par(mfrow = c(1,3))
plot(globe1, main = 'm9'); title(main = 'globe1')
plot(globe2, main = 'm9'); title(main = 'globe2')
plot(globe3, main = 'm9'); title(main = 'globe3')



## --------------------------------------------------------------
set.seed(2)
mesh = inla.mesh.2d(x_loc, max.edge = c(0.1,0.5))
spde = inla.spde2.matern(mesh, loc = x_loc)
Q = inla.spde.precision(spde, theta = c(0,0))
sample = inla.qsample(n = 2, Q)


## --------------------------------------------------------------
proj <- inla.mesh.projector(mesh, dims = c(100, 100))
proj2 <- inla.mesh.projector(mesh, dims = c(50, 50))
proj3 <- inla.mesh.projector(mesh, dims = c(10, 10))


## --------------------------------------------------------------
sample_proj = inla.mesh.project(proj, field = sample[,1])
sample_proj2 = inla.mesh.project(proj2, field = sample[,1])
sample_proj3 = inla.mesh.project(proj3, field = sample[,1])


## --------------------------------------------------------------
par(mfrow = c(1,3))
image(proj$x, proj$y, sample_proj , xlim = c(0,1), ylim = c(0,1),
      xlab = '',ylab = '')
contour(proj$x, proj$y, sample_proj, add = T)

image(proj2$x, proj2$y, sample_proj2, xlim = c(0,1), ylim = c(0,1),
      xlab = '',ylab = '')
contour(proj2$x, proj2$y, sample_proj2, add = T)

image(proj3$x, proj3$y, sample_proj3, xlim = c(0,1), ylim = c(0,1),
      xlab = '',ylab = '')
contour(proj3$x, proj3$y, sample_proj3, add = T)



## --------------------------------------------------------------
par(mfrow = c(1,2))
set.seed(123)
mesh2 <- inla.mesh.create(globe = 10)
spde.glob = inla.spde2.matern(mesh2, loc = mesh2$loc)
Q.glob = inla.spde.precision(spde.glob, theta = c(0, 0))
x.glob = inla.qsample(n = 2, Q.glob)

proj2a <- inla.mesh.projector(mesh2, projection = "longlat",
                                   dims = c(361, 181))
proj2b <- inla.mesh.projector(mesh2, projection = "mollweide",
                                   dims = c(361, 181))

image(proj2b$x, proj2b$y, inla.mesh.project(proj2b, field = x.glob[,1]), main = 'Mollwide Projection', xlab = '', ylab = '')
contour(proj2b$x, proj2b$y, inla.mesh.project(proj2b, field = x.glob[,1]), add = T)

image(proj2a$x, proj2a$y, inla.mesh.project(proj2a, field = x.glob[,1]), main = 'Long-Lat projection', xlab = '', ylab = '')
contour(proj2a$x, proj2a$y, inla.mesh.project(proj2a, field = x.glob[,1]), add = T)


