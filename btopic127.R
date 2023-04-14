## ----setup, include=FALSE-----------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- message = F, warning=FALSE----------------------------------
library(INLA)
library(rgdal)
library(rgeos)
library(ggmap)
set.seed(2018)


## -----------------------------------------------------------------
# First polygon
pl1 <- Polygon(cbind(c(0,15,15,0,0), c(5,0,20,20,5)), hole=FALSE)
# Hole in the first polygon
h1 <- Polygon(cbind(c(5,12,10,5,5), c(7,7,15,15,7)), hole=TRUE)
# Second polygon
pl2 <- Polygon(cbind(c(15,20,20,30,30,15,15), c(10,10,0,0,20,20,10)), hole=FALSE)


## -----------------------------------------------------------------
sp <- SpatialPolygons(list(Polygons(list(pl1, h1), '0'), 
                           Polygons(list(pl2), '1'))) 


## -----------------------------------------------------------------
plot(sp, main = 'Polygon')


## -----------------------------------------------------------------
mesh <- inla.mesh.2d(boundary=sp, max.edge=2)
plot(mesh, main = 'Mesh')


## -----------------------------------------------------------------
mesh2 <- inla.mesh.2d(boundary=sp, max.edge=c(1,4))


## ---- echo = F----------------------------------------------------
plot(mesh2, main = 'Extended Mesh')
title(main = 'Extended Mesh')


## ---- message=FALSE,warning=FALSE,eval=FALSE----------------------
## shape <- readOGR(dsn = ".", layer = "simplified_land_polygons")


## ----eval=FALSE---------------------------------------------------
## png("fig/btopic127/shapeplot1.png", width = 480*4, height=480*4)
## plot(shape)
## dev.off()
## 
## # Open the saved file to see the following plot


## ---- eval = F----------------------------------------------------
## shape2 <- spTransform(shape, CRS("+proj=longlat +datum=WGS84"))


## ---- eval = F----------------------------------------------------
## png("fig/btopic127/shapeplot2.png", width = 480*4, height=480*4)
## plot(shape2)
## dev.off()


## ---- eval = F----------------------------------------------------
## png("fig/btopic127/shapeplot3.png", width = 480*4, height=480*4)
## plot(shape2, xlim=c(22.80, 26.85), ylim=c(34.88, 38.06))
## dev.off()


## -----------------------------------------------------------------
# You can run this code chunk
pl1 <- Polygon(cbind(c(22.80, 22.80, 26.85, 26.85, 22.80), 
                     c(38.06, 34.88, 34.88, 38.06, 38.06)), hole=FALSE)
sp <- SpatialPolygons(list(Polygons(list(pl1), '0')), 
                      proj4string =CRS("+proj=longlat +datum=WGS84"))


## ---- eval = F----------------------------------------------------
## shape3 <- gBuffer(shape2, byid=TRUE, width=0) #for gIntersection operation
## shape4 <- gIntersection(shape3, sp) #this is the heavy guy
## shape5 = gSimplify(shape4, tol=0.001) #use tol to control the precision


## ---- eval = F----------------------------------------------------
## shape_df = as(shape5, 'SpatialPolygonsDataFrame')
## writeOGR(shape_df, layer = 'shape5','temp/' ,driver="ESRI Shapefile")


## -----------------------------------------------------------------
shape5 <- readOGR(dsn = "data/btopic127/", layer = "shape5")


## ---- message = F, fig.width=15, fig.height=15--------------------
# create a polygon representing the area of interest
pll = Polygon(sp@polygons[[1]]@Polygons[[1]]@coords, hole = F)

# number of polygons representing the islands
n_poly = length(shape5@polygons[[1]]@Polygons)
idx = seq(1:n_poly)

# create a list of holes
hole_list = lapply(idx, function(n) Polygon(shape5@polygons[[1]]@Polygons[[n]]@coords, hole = T)) 

# create the final Spatial Polygons
new_sp = SpatialPolygons(list(Polygons(append(hole_list, pll),'1')))

# take a look
plot(new_sp)


## ---- message = F, fig.width=8, fig.height=8----------------------
mesh_poly = inla.mesh.2d(boundary = new_sp, max.edge = 0.2)
plot(mesh_poly, main = '')


## ---- message = F, fig.width=8, fig.height=8----------------------
mesh_poly2 = inla.mesh.2d(boundary = new_sp, max.edge = c(0.1, 0.4),
                          cutoff = c(0.02))
plot(mesh_poly2, main = '')


## ---- message = F, fig.width=6, fig.height=6, warning=F-----------
try({
  glgmap = get_map(location = shape5@bbox, maptype= "terrain", zoom = 7)
  p = ggmap(glgmap) #p is the plot of glgmap
  print(p)
})


## ---- message = F, fig.width=6, fig.height=6----------------------
try({
  p1 = p + xlab("Longitude") + ylab("Latitude") 
  p2 = p1 + geom_polygon(data = fortify(shape5),
                       aes(long, lat, group = group),
                       fill = "orange", colour = "red", alpha = 0.2)
  p3 = p2 + geom_polygon(data = fortify(sp),
                       aes(long, lat, group = group),
                       fill = NA, colour = "black", alpha = 1)
  print(p3)
})



## -----------------------------------------------------------------
#```{r, message = F, fig.width=8, fig.height=8}
## 13th april 2023:
## This code no longer runs! Bummer!
## If anyone knows how to fix, please let me know!
try({
  shape_utm = spTransform(shape5, CRS("+proj=utm +zone=35S ellps=WGS84")) 
  plot(shape_utm)
})

