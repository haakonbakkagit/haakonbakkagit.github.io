### Model 121 as in the theorem in the manuscript 

### space-time model 121
stmodel121.interpret <- function(cmd=c("graph", "Q", "mu", "initial", "log.norm.const", 
                             "log.prior", "quit"), theta = NULL, args = NULL) {

  stopifnot(length(obj$lambdas) == 3)
  stopifnot(inherits(obj$mesh.space, "inla.mesh"))
  stopifnot(inherits(obj$mesh.time, "inla.mesh.1d"))
  
  ## Set fixed.theta if it exists
  fixed.theta = NULL
  fixed.theta = obj$fixed.theta
  
  ## Spatial and temporal FEM matrices
  sfe = inla.mesh.fem(obj$mesh.space, order = 4)
  tfe = inla.mesh.fem(obj$mesh.time, order = 2)
  
  ## Temporal matrices according to our notation
  M0 = tfe$c0
  stopifnot(abs(M0[1,1] - 0.5*M0[2,2])<1e-3)
  N = nrow(M0)
  M1 = sparseMatrix(i=c(1,N), j=c(1,N), x=0.5)
  M2 = tfe$g1
  
  interpret.theta <- function(n, theta) {
    if (is.null(fixed.theta)) {
      ## Assume the input theta is log(rt, rs, sigma)
      theta.interpret = theta
    } else {
      stopifnot(length(fixed.theta)<4)
      theta.interpret = fixed.theta
      theta.interpret[is.na(fixed.theta)] = theta
    }
    ## Note that theta is log(range_t, range_s, sigma)
    alpha_t = 1; alpha_s = 2; alpha_e = 1;
    alpha = alpha_e + alpha_s*(alpha_t-1/2);
    nu.s = alpha -1; nu.t = alpha_t - 0.5;
    ## old&davids: 
    #c1 = gamma(alpha_t - 1/2)*gamma(alpha-1)/(gamma(alpha_t) *gamma(alpha) *4*sqrt(pi))
    ## New
    c1 = gamma(alpha_t - 1/2)*gamma(alpha-1)/(gamma(alpha_t) *gamma(alpha) *8*pi^1.5)
    
    ## Define the log-gamma's
    theta.gam = rep(NA, 3)
    theta.gam[2] = 0.5*log(8*nu.s) - theta.interpret[2]
    theta.gam[1] = theta.interpret[1] - 0.5*log(8*(alpha_t-1/2)) + alpha_s * theta.gam[2]
    theta.gam[3] = 0.5*log(c1) - 0.5*theta.gam[1] - (alpha-1)*theta.gam[2] - theta.interpret[3]
    return(theta.gam)
  }
  graph <- function(n, theta) {
    require(methods)
    ntheta = length(initial())
    G1 = Q(theta = (1:ntheta)/3.217233456)
    G1[G1 != 0] = 1
    G2 = Q(theta = (1:ntheta)^2/12.1543534)
    G2[G2 != 0] = 1
    return(G1 + G2)
  }
  
  
  Q <- function(n, theta) {
    theta.gam <- interpret.theta(theta=theta)
    gt <- exp(theta.gam[1]) ## squared \gamma_t
    gs <- exp(theta.gam[2]) ## \gamma_s
    ge2 <- exp(2*theta.gam[3]) ## squared \gamma_e
##    cat("log: gt, gs, ge", theta.gam, '\n') 
    
    
    Q = (kronecker(gt^2*M2,
                   gs^2*sfe$c0 + sfe$g1) +
           kronecker(M0,
                     gs^6*sfe$c0 +
                       gs^4*sfe$g1 +
                       gs^2*sfe$g2 +
                       sfe$g3) +
           kronecker(2*gt*M1,
                     gs^4*sfe$c0 +
                       2*gs^2*sfe$g1 +
                       sfe$g2 )) * ge2
    return(Q)
  }
  
  mu <- function(n, theta) return(numeric(0))
  log.norm.const <- function(n, theta) return(numeric(0))
  ## pc priors with 3 lambdas
  log.prior <- function(n, theta) { 
    if (is.null(fixed.theta)) {
      ## Assume the input theta is log(rt, rs, sigma)
      #theta = theta
    } else {
      ## In this case the prior will not be properly re-scaled
      ## But that makes no difference to INLA (except mlik etc.)
      stopifnot(length(fixed.theta)<4)
      theta.interpret = fixed.theta
      theta.interpret[is.na(fixed.theta)] = theta
      theta = theta.interpret
    }
    ## lambdas for the corresponding thetas are in obj
    lambdas = obj$lambdas
    
    ## log prior value
    #i = 1 # log(rt) # pc prior for range in d=1
    val = log(lambdas[1]) - lambdas[1] * exp(-0.5*theta[1]) + 
      log(0.5) - 0.5*theta[1]
    #i = 2 # log(rs) # pc prior for range in d=2
    val = val + log(lambdas[2]) - lambdas[2] * exp(-theta[2]) + 
      -theta[2]
    #i = 3 # sigma # pc prior for sigma
    val = val + log(lambdas[3]) - lambdas[3] * exp(theta[3]) + 
      theta[3]
    
    return(val)
  }
  initial <- function(n, theta) {
    if (is.null(fixed.theta)) {
      return(c(0,0,0))
    } else {
      return(rep(0, sum(is.na(fixed.theta))))
    }
  }
      
  quit <- function(n, theta) return(invisible())
  cmd <- match.arg(cmd)
  val <- do.call(cmd, args=list(n=as.integer(args$n), theta=theta))
  return(val)
}


