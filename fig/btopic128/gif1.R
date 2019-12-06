library(INLA)
library(fields)

## ----build the shape of the barrier-----------------------------

set.seed(2017)
set.inla.seed = 2017

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

## ----build the mesh----------------------------------------------

max.edge.length = 0.2#改小一点

n=4
B=c(0,0, 10,0, 0,10, 10,10)

for (i in seq(3, 7, 0.3)){
  
  B=c(B,5,i)
  n=n+1
  
}

loc1 = matrix(B, n, 2, byrow = T)#把要画点的地方加进去

seg = inla.sp2segment(poly.original)

mesh = inla.mesh.2d(loc=loc1, interior = seg, 
                    max.e = max.edge.length, offset=1)

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

plot(mesh, main="Mesh and Omega")
plot(poly.barrier, add=T, col='lightblue')
plot(mesh, add=T)
points(loc1)

local.plot.field = function(field, ...){
  xlim = c(2, 8); ylim = xlim;
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 grid 
  #   for plots
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, ...)  
  # - Use image.plot to get nice colors and legend
}
print(mesh$n)

## ----define the barrier model--------------------------------------

barrier.model = inla.barrier.pcmatern(mesh, barrier.triangles = barrier, prior.range = c(1.44, 0.5), prior.sigma = c(0.7, 0.5), range.fraction = 0.1)

range = 1.2
# - the spatial range parameter
Q = inla.rgeneric.q(barrier.model, "Q", theta = c(0, log(range)))

u = inla.qsample(n=1, Q=Q, seed = set.inla.seed)
u = u[ ,1]
# - access the first sample

local.plot.field(u, main="The true (simulated) spatial field")
plot(poly.barrier, add=T, col='grey')

##----calculate the prior distribution-------------------------------

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

corr = local.find.correlation(Q, loc = c(5,3), mesh)

local.plot.field(corr, mesh, zlim=c(0.1, 1))
plot(poly.barrier, add=T, col="grey")

## ----create and save the GIF---------------------------------------

library("animation")
saveGIF(
  {
    
    
    oopt = ani.options(interval = 0.05, nmax =100)
    
    #create the animation
    
    for (i in seq(3, 7, 0.3)){
      
      corr = local.find.correlation(Q, loc = c(5,i), mesh)
      
      local.plot.field(corr, mesh, zlim=c(0.1, 1))
      plot(poly.barrier, add=T, col="grey")
      
      #wait the time period set by interval
      
      ani.pause()
      
    }
    
    #reload the animation options
    
    ani.options(oopt)
    
  },movie.name="prior of barrier model1.gif",img.name="prior1")

