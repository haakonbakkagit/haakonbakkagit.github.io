### DESCRIPTION ###
# These functions define the Barrier and Different Terrains models


### For internal use ###
# This file should not be edited!
# - Only perfect backwards compatible changes allowed
# - Last update: 5. April 2017
# A new filename should be created, and all files with the current filename should be systematically edited and verified.
# (Master-Slave-document organisation was abandoned)

### LOAD LIBRARIES ###
library(INLA)
library(rgeos)
library(fields)
library(rgdal)

### GENERAL FUNCTIONS ###
bakka.square.polygon = function(xlim, ylim, ret.SP = F, id=runif(1)){
  # - ret.SP=T : Return a SpatialPolygons object
  # - ret.SP=F : Return a Polygon object
  xlim = range(xlim); ylim = range(ylim)
  corner1 = c(xlim[1], ylim[2])
  corner2 = c(xlim[2], ylim[1])
  poly = Polygon(rbind(corner1, c(corner1[1], corner2[2]), corner2, c(corner2[1], corner1[2]), corner1), hole = FALSE)
  if (ret.SP) {
    return(SpatialPolygons(list(Polygons(list(poly), ID = id))))
  } else {
    return(poly)
  }
}

dt.mesh.addon.posTri <- function(mesh){
  # - Add on two attributes to the mesh object
  mesh$t = length(mesh$graph$tv[,1])
  mesh$posTri = matrix(0, mesh$t, 2)
  for (t in 1:mesh$t){
    temp = mesh$loc[mesh$graph$tv[t, ], ]
    mesh$posTri[t,] = colMeans(temp)[c(1,2)] 
  }
  return(mesh)
}

dt.Omega <- function(list_of_subdomains, mesh){
  # - list_of_subdomains must be a list of numeric 'c(...)'
  # - This fun creates a legal Omega from just a list of subdomains
  # - It does not make sure every triangle is included somewhere
  # - But it does make sure there are no multiplicities
  
  if (length(list_of_subdomains[[1]])==0){
    stop('No proper first subdomain')
  }
  Omega = list()
  Omega[[1]] = list_of_subdomains[[1]]
  if (length(list_of_subdomains)==1){
    remainder = setdiff(1:mesh$t, Omega[[1]])
    if (length(remainder)>0){
      Omega[[length(Omega)+1]] = remainder
    }
    return (Omega)
  }
  usedTriangles = Omega[[1]]
  for (i in 2:length(list_of_subdomains)){
    Omega[[i]] = setdiff(list_of_subdomains[[i]], usedTriangles)
    usedTriangles = union(usedTriangles, Omega[[i]])
  }
  ## Remove empty Omega parts
  for (i in 1:length(list_of_subdomains)){
    if (length(Omega[[length(Omega)]])==0) Omega[[length(Omega)]] = NULL
  }
  
  ## Add on remainder
  remainder = setdiff(1:mesh$t, usedTriangles)
  if (length(remainder)>0){
    Omega[[length(Omega)+1]] = remainder
  }
  return (Omega)
}


dt.polygon.omega = function (mesh, Omega) {
  # - constructs SpatialPolygons for the different subdomains (areas)
  stopifnot(class(mesh) == 'inla.mesh')
  # - requires an inla mesh to work
  
  Omega.SP.list = list()
  for (j in 1:length(Omega)) {
    poly.list = list()
    for (tri in Omega[[j]]){
      px = mesh$graph$tv[tri, ]
      temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
      poly.list = c(poly.list , Polygon(rbind(temp[ ,1:2], temp[1, 1:2]), hole=F))
    }
    mesh.polys = SpatialPolygons(list(Polygons(poly.list, ID='noid')))
    Omega.SP.list[[j]] = gUnaryUnion(mesh.polys)
  }
  return(Omega.SP.list)
}


dt.nearest.mesh.index = function (location, mesh, verbose = T) {
  # Uses the A-construction and picks one of the nearest mesh nodes
  
  A.tmp = inla.spde.make.A(mesh=mesh, loc = matrix(c(location[1], location[2]),1,2))
  node.id = which.max(A.tmp[1, ])

  return(node.id)
}



### PRECISION MATRIX FUNCTIONS ###
# I.e. Solve the differential equation (the SPDE)

dt.precision.new <- function(spde, ranges, sigma=1){
  # - This function computes a specific precision matrix
  # - spde is the Different Terrains model
  # - spde contains all the needed matrices to solve the SPDE
  # - the ranges and sigma are the hyperparameters that determine Q
  
  if(is.null(ranges)) stop("ranges cannot be NULL")
  if(any(is.na(ranges))) stop("No range can be NA")
  
  xi = length(ranges)
  if (xi != length(spde$D)){
    print('dt.precision has encountered an error. Will stop.')
    stop ('Ranges do no correspond to spde')
  }
  if (any(ranges < 0.001)){
    warning('This hyper parameter value will probably fail. A very small maximum edge length needed in the mesh.')
  }
  
  Cdiag = ranges[1]^2* spde$C[[1]] # already raised to power 2
  if (xi > 1){
    for (k in 2:xi){
      Cdiag = Cdiag + ranges[k]^2*spde$C[[k]]
    }
  }
  N = length(Cdiag)
  Cinv = sparseMatrix(i=1:N, j = 1:N, x=1/Cdiag, dims = c(N,N), giveCsparse=FALSE) 
  
  A = spde$I  
  for (k in 1:xi){
    A = A+ (ranges[k]^2/8)*spde$D[[k]]
  }
  
  Q = t(A)%*%Cinv%*%A*(1/sigma^2)/(pi/2) *4 
  Q = inla.as.dgTMatrix(Q)
  return (Q)
}


dt.fem.matrices <- function(mesh, Omega){
  # - This function computes the Finite Element matrices
  # - - this is needed to compute the precision matrix Q later
  
  xi = length(Omega)
  spde = list()
  spde$I = dt.fem.identity(mesh)
  spde$D = list()
  spde$C = list()
  for (k in 1:xi) {
    spde$D[[k]] = dt.fem.laplace(mesh, Omega[[k]])
  }
  for (k in 1:xi){
    spde$C[[k]] = dt.fem.white(mesh, Omega[[k]])
  }
  spde$hdim = xi
  return(spde)
}


dt.fem.white <- function(mesh, subdomain){
  # - This function computes the Finite Element matrix of the white noise on ONE subdomain (area)
  # - This matrix is a diagonal matrix
  
  ## Pre-allocation
  Ck = rep(0, mesh$n)

  for (t in subdomain){
    ## Node indexes for triangle t:
    px = mesh$graph$tv[t, ]
    temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
    p1 = t(t(temp[1, c(1,2)]))
    p2 = t(t(temp[2, c(1,2)]))
    p3 = t(t(temp[3, c(1,2)]))
    Ts = cbind(p2-p1, p3-p1) # is the transformation to reference triangle
    area=  abs(det(Ts)) * 0.5
    for (i in 1:3){
      Ck[px[i]] = Ck[px[i]] + area
    }
  }
  return (Ck)
}


dt.fem.identity <- function(mesh){
  # - this function computes the Finite Element matrix for the '1' in the SPDE (the identity operator)
  # - this operator does not depend on the subdomains
  
  ## Preallocation
  len = length(mesh$graph$tv[,1])
  index.i = rep(0,len * 6)
  index.j = rep(0,len * 6)
  Aij = rep(0,len * 6) # matrix values
  counter = 1
  
  for (t in 1:len){
    ## Node indexes for triangle t:
    px = mesh$graph$tv[t, ]
    temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
    p1 = t(t(temp[1, c(1,2)]))
    p2 = t(t(temp[2, c(1,2)]))
    p3 = t(t(temp[3, c(1,2)]))
    Ts = cbind(p2-p1, p3-p1) # is the transformation to reference triangle
    twiceArea=  abs(det(Ts))
    
    for (i in 1:3){
      index.i[counter] = px[i]
      index.j[counter] = px[i] # same
      Aij[counter] = (twiceArea)*1/12
      counter = counter + 1
    }
    
    for (i in 1:2){
      for (j in (i+1):3)
        index.i[counter] = px[i]
      index.j[counter] = px[j]
      Aij[counter]= (twiceArea)*1/24
      counter=counter + 1
      # symmetry:
      index.i[counter] = px[j]
      index.j[counter] = px[i]
      Aij[counter]= (twiceArea)*1/24
      counter=counter + 1
    }
  }
  
  I = sparseMatrix(i=index.i, j = index.j, x=Aij, dims = c(mesh$n, mesh$n), giveCsparse=FALSE) 

  return (I)
}



dt.fem.laplace <- function(mesh, subdomain){
  # - This function computes the Finite Element matrix of the laplace operator on ONE subdomain (area)
  # - This matrix is very sparse
  
  # The nabla phi's
  Nphix = rbind(c(-1,-1), c(1,0), c(0,1))
  len = length(subdomain)
  index.i = rep(0,len * 9)
  index.j = rep(0,len * 9)
  Aij = rep(0,len * 9) # matrix values
  counter = 1
  
  for (tri in subdomain){
    px = mesh$graph$tv[tri, ]
    temp = mesh$loc[px, ] # is a 3 by 3 matrix of node locations
    p1 = t(t(temp[1, c(1,2)]))
    p2 = t(t(temp[2, c(1,2)]))
    p3 = t(t(temp[3, c(1,2)]))
    Ts = cbind(p2-p1, p3-p1) # is the transformation to reference triangle
    TTTinv = solve(t(Ts)%*%Ts)
    area=  abs(det(Ts)) * 0.5
    
    for (k in 1:3){
      for (m in 1:3){
        tmp = (3*m+k-4)*length(subdomain)
        index.i[(tmp + counter)] = px[k]
        index.j[(tmp + counter)] = px[m]
        Aij[(tmp + counter)] = area*Nphix[k, c(1,2) ]%*%TTTinv%*%as.matrix(Nphix[m, c(1,2) ])
      }
    }
    counter = counter + 1
  }
  
  Dk = sparseMatrix(i=index.i, j = index.j, x=Aij, dims = c(mesh$n, mesh$n), giveCsparse=FALSE) 
  return (Dk)
}



### FUNCTIONS FOR RGENERIC MODEL SETUP ###

## Check
if (!exists('dt.precision.new')) {
  stop('Error: You need additional functions')
}

dt.create.Q = function(mesh, Omega, initial.theta = NULL, 
                       fixed.ranges = NULL, copy.ranges = NULL, copy.ranges.frac = NULL)
{
  # This function creates the precision matrix with environment
  
  ## Input
  # fixed.ranges
  # - NULL, or a list of values for the fixed and NA (for the rest)
  # initial.theta
  # - Initial values for setting the hypers 
  # - Internal scale: log(sigma) log(ranges)
  # - will be overridden by control.mode in the inla(...) statement
  
  ## Output Q
  # - is a precision function with environment
  # - contains n, ntheta, spde, dt.precision.new, initial.theta, and possibly others
  
  ## Now implemented: copy.ranges
  # - NULL or c(values) list of values smaller than index for an index to copy
  # - The purpose is to be able to use this to run a stationary model without creating another spde variable
  # - copy.ranges.frac gives the fraction to multiply the range with
  
  ## Copy variables to current environment/scope
  dt.precision.new = dt.precision.new
  initial.theta = initial.theta
  fixed.ranges = fixed.ranges
  copy.ranges = copy.ranges
  copy.ranges.frac = copy.ranges.frac
  
  ## Create spde object (time consuming)
  spde = dt.fem.matrices(mesh = mesh, Omega = Omega)
  
  ## Create Q function
  if (is.null(fixed.ranges)){
    
    if (is.null(copy.ranges)) {
      Q = function(theta){
        return(dt.precision.new(spde=spde, ranges=exp(theta[-1]), sigma=exp(theta[1])))
      }
      
      ntheta = 1+length(spde$D)
      # - create variable needed for rgeneric
      # - ntheta means number of thetas  
    } else {
      ## Check input and use default values
      if(length(copy.ranges)!=length(spde$D)) stop("Wrong input length for copy.ranges")
      if (is.null(copy.ranges.frac)) copy.ranges.frac=rep(1, length(spde$D))
      if(length(copy.ranges.frac)!=length(spde$D)) stop("Wrong input length for copy.ranges.frac")
      copy.ranges.frac[is.na(copy.ranges.frac)] = 1
      if(any(copy.ranges.frac[is.na(copy.ranges)] != 1)) warning("copy.ranges.frac ignored where you are not copying any ranges")
      
      if(!all(is.na(copy.ranges[copy.ranges[!is.na(copy.ranges)]]))) stop("You can only copy ranges that themselves are NA")
      
      ## Construct the Q function
      Q = function(theta){
        id.na = which(is.na(copy.ranges))
        id.notna = which(!is.na(copy.ranges))
        
        exptheta = exp(theta[-1])
        ranges = rep(NA, length(copy.ranges))
        ranges[id.na] = exptheta 
        # - this dimension should be correct
        ranges[id.notna] = ranges[copy.ranges[id.notna]]*copy.ranges.frac[id.notna]
        # - this is copying from the other ranges values instead of from exptheta
        # - only the appropriate scalings are taken into account
        
        return(dt.precision.new(spde=spde, ranges=ranges, sigma=exp(theta[1])))
      }
      
      ntheta = 1+sum(1*is.na(copy.ranges))
      
    }
    
  } else {
    if (!is.null(copy.ranges)) stop('Cannot have both fixed and copied ranges')
    
    # we have to fix some hyperparameters
    stopifnot(length(fixed.ranges) == length(spde$D))
    Q = function(theta){
      if (length(theta) != sum(is.na(fixed.ranges))+1){
        print(theta); print(fixed.ranges); stop("Num ranges not equal to num not fixed ranges")
      }
      internal.ranges = fixed.ranges
      internal.ranges[is.na(internal.ranges)] = exp(theta[-1])
      return(dt.precision.new(spde, internal.ranges, sigma = exp(theta[1])))
    }
    ntheta = 1+sum(is.na(fixed.ranges))
    # - create variable needed for rgeneric
    # - ntheta means number of thetas
  }
  
  n = dim(spde$I)[1]

  if(is.null(initial.theta)) {
    initial.theta = rep(0, ntheta)
  }

  return(Q)
  
}

dt.validate.model.component = function(Q) {
  # Q is a function with an environment
  e = environment(Q)
  if (e$ntheta != length(e$initial.theta)) stop("wrong ntheta or initial")
  
  stopifnot(any(names(e)=="n"))
  stopifnot(any(names(e)=="ntheta"))
  stopifnot(any(names(e)=="spde"))
  stopifnot(any(names(e)=="dt.precision.new"))
  stopifnot(any(names(e)=="initial.theta"))
  
  tmp = Q(e$initial.theta)
  stopifnot(all(dim(tmp==e$n)))
  
  return(TRUE)
}

dt.create.prior.log.exp = function (prior.param) {
  # This is the log-prior for the internal parametrisation theta
  # Both log of probability and log of exponential dist
  # theta = log(sigma), log(range1), log(range2), ...
  
  ## Input
  # parameters are the lambdas in the exponential distribution
  # - the first is for sigma
  # - the second is for all the ranges
  
  ## Move to current scope (environment)
  prior.param = prior.param
  
  log.prior = function(theta) {
    lambda0 = prior.param[1]
    lambda1 = prior.param[2]
    ntheta = length(theta)
    
    ## Prior for standard deviation
    val = 0 + log(lambda0) - lambda0*exp(theta[1]) + theta[1]
    
    ## Prior for range(s)
    for (i in 2:ntheta) {
      val = val + log(lambda1) - lambda1*exp(-theta[i]) + -theta[i]
    }
    return(val)
  }
  
  return(log.prior)
  # - this environment includes the prior parameters
}


dt.inla.model = function(Q, log.prior) {

  ## Input
  # Q and log.prior
  # - must be functions of theta
  # - must be self-contained, with constants fixed in environment
  
  stopifnot(dt.validate.model.component(Q))
  
  ## Input verification
  if(identical(environment(Q),.GlobalEnv)) stop("Cannot use global environments")
  if(identical(environment(log.prior),.GlobalEnv)) stop("Cannot use global environments")
  if(anyDuplicated(ls(environment(Q)), ls(environment(log.prior)))) stop("Duplicate variables")
  
  model.rgeneric = inla.rgeneric.define(model = general.rgeneric.model, Q=Q, log.prior=log.prior)
  
  return(model.rgeneric)
}


general.rgeneric.model = function(
  # - this function is the model component definition in the rgeneric inla framework
  # - see demo(rgeneric)
  # - see inla.doc('rgeneric')
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
          "log.prior", "quit"),
  theta = NULL)
  
{
  ## Input
  # theta : only used when constructing Q
  
  ## Assumed functions
  # Q(theta)
  # log.prior(theta)
  # These are already defined
  
  ## New functions
  graph = function(theta) 
  {
    require(methods)
    ntheta = length(initial())
    
    G1 = Q(theta=(1:ntheta)/3.217233456)
    G1[G1 != 0] = 1
    G2 = Q(theta=(1:ntheta)^2/12.1543534)
    G2[G2 != 0] = 1
    
    return (G1+G2)
    # - this should never give zeroes by random.
  }
  
  
  log.norm.const = function(theta) numeric(0)
  
  initial = function(theta)
  {
    ## Switch:
    #print.v = print
    print.v = function(x) { }
    
    print.v("Test hi - this should work!")
    print.v(ls())
    print.v("next one up 1")
    print.v(ls(envir = (parent.env(environment()))))
    print.v("next one up 2")
    print.v(parent.env(parent.env(environment())))
    print.v(ls(envir = parent.env(parent.env(environment()))))
    print.v("Q environment")
    print.v(ls(environment(Q)))
    
    return (environment(Q)$initial.theta) 
  }
  
  quit <- function(theta) invisible()
  
  mu <- function(theta) numeric(0)
  
  val = do.call(match.arg(cmd), args = list(theta))
  return (val) 
} # end function definition






