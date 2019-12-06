## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE---------------------------------------
library(INLA)
rm(list=ls())
options(width=70, digits=2)
set.seed(2017)


## ------------------------------------------------------------------------
loc = matrix(c(1, 0, 1, 1, 0, 1, 0,0), ncol=2, byrow = T)
mesh1 = inla.mesh.2d(loc, max.edge = 0.2)
plot(mesh1)


## ------------------------------------------------------------------------
loc = matrix(c(1, 0, 1, 1, 0, 1, 0,0), ncol=2, byrow = T)
mesh2 = inla.mesh.2d(loc, max.edge = 0.2, offset = c(0, 0.05))
plot(mesh2)


## ------------------------------------------------------------------------
loc = matrix(c(1, 0, 1, 1, 0, 1, 0,0), ncol=2, byrow = T)
mesh3 = inla.mesh.2d(loc, max.edge = c(0.2, 0.2), offset = c(0, 0.5))
plot(mesh3)


## ------------------------------------------------------------------------
loc = matrix(c(1, 0, 1, 1, 0, 1, 0, 0, 0.5, 0.5, 0.499, 0.499, 0.3, 0.3, 0.2999, 0.2999), ncol=2, byrow = T)
mesh4 = inla.mesh.2d(loc, max.edge = 0.2, cutoff = 1e-12)
plot(mesh4)
points(loc)

