## ----setup, include=FALSE-----------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE--------------------------------
library(INLA); library(sp)

set.seed(2016)
set.inla.seed = 2016


## -----------------------------------------------------------------
## Load data
load(file = "data/WebSiteData-Archipelago.RData")
# - if you have saved the file locally

## What is loaded
# - poly.water is our study area
# - df is our dataframe to be analysed
# - dat is the orginial dataframe
str(poly.water, 1)


## -----------------------------------------------------------------
summary(poly.water)
plot(poly.water, axes=T)
points(df$locx, df$locy, col="blue", cex=0.5)


## -----------------------------------------------------------------
summary(df)


## -----------------------------------------------------------------
max.edge = 0.95
# - some chosen constant
# - results should not be sensitive to this (if you have a good mesh)
# - max.edge = diff(range(df$locx))/15
mesh1 = inla.mesh.2d(loc=cbind(df$locx, df$locy),
                    max.edge = max.edge)
plot(mesh1, main="1st attempt"); points(df$locx, df$locy, col="blue")


## -----------------------------------------------------------------
max.edge = 0.95
# - as before
bound.outer = 4.6
# - the outer boundary I want to use for my mesh
# - some chosen constant
# - results should not be sensitive to this
# - bound.outer = diff(range(df$locx))/3
mesh2 = inla.mesh.2d(loc=cbind(df$locx, df$locy),
                    max.edge = c(1,5)*max.edge,
# - use 5 times max.edge in the outer extension/offset/boundary
                    cutoff = max.edge/5,
                    offset = c(max.edge, bound.outer))
plot(mesh2, main="2nd attempt"); points(df$locx, df$locy, col="blue")


## ---- dpi=72*1.5--------------------------------------------------
max.edge = 0.95
# - as before
bound.outer = 4.6
# - as before
mesh3 = inla.mesh.2d(boundary = poly.water,
                     loc=cbind(df$locx, df$locy),
                    max.edge = c(1,5)*max.edge,
# - use 5 times max.edge in the outer extension/offset/boundary
                    cutoff = max.edge/5,
                    offset = c(max.edge, bound.outer))
plot(mesh3, main="3rd attempt"); points(df$locx, df$locy, col="red")


## ---- include=FALSE-----------------------------------------------
## Backup plan to increase ppi
#ppi = 300; png("mygraph.png", )
#dev.args=list(width=6*ppi, height=6*ppi, res=ppi)


## ---- dpi=72*2----------------------------------------------------
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


## -----------------------------------------------------------------
mesh4$n


## -----------------------------------------------------------------
in.water = over(poly.water, SpatialPoints(cbind(df$locx, df$locy)), returnList=T)[[1]]
print(paste("There are", nrow(df)-length(in.water), "points on land in the original polygon"))

