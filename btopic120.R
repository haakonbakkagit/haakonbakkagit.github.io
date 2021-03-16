## ----setup, include=FALSE, warning=FALSE-----------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- message=FALSE--------------------------------------------
library(INLA)


## --------------------------------------------------------------
N = 10
rho = 0.95
Q = matrix(0, N, N)
diag(Q) = 1+rho^2
for (i in 1:(N-1)) {
  Q[i, i+1] = -rho
  Q[i+1, i] = -rho
}
print(Q)


## --------------------------------------------------------------
Q[1,1] = 1
Q[N,N] = 1


## --------------------------------------------------------------
precision.ar1 = function(N, rho){
  Q = matrix(0, N, N)
  diag(Q) = 1+rho^2
  for (i in 1:(N-1)) {
    Q[i, i+1] = -rho
    Q[i+1, i] = -rho
  }
  Q[1,1] = 1
  Q[N,N] = 1
  return(Q)
}


## --------------------------------------------------------------
Q = precision.ar1(10, 0.9)
print(Q)


## --------------------------------------------------------------
Q = precision.ar1(10, 0.9)
Q = as(Q, "sparseMatrix")
print(Q)


## --------------------------------------------------------------
str(Q)


## --------------------------------------------------------------
for (N in c(10, 1E2, 1E3, 5E3 )) {
  Q = precision.ar1(N, 0.9)
  os1 = round(object.size(Q)/1000)
  Q = as(Q, "sparseMatrix")
  os2 = round(object.size(Q)/1000)
  print(paste0("For N is ", N, " we go from ", os1, " kb to ", os2, " kb"))
}


## --------------------------------------------------------------
Q = precision.ar1(1000, 0.99)
L = chol(Q)
print(sum(abs(Q- t(L)%*%L)))


## --------------------------------------------------------------
pt1 = rep(0.1, nrow(Q))
log.prob.pt1 = -nrow(Q)/2*log(2*pi) + 0.5*2*sum(log(diag(L))) - 0.5 *pt1 %*% Q %*% pt1
print(log.prob.pt1)


## --------------------------------------------------------------
set.seed(2017)
z.sample = rnorm(nrow(Q))
u.sample = solve(L, z.sample)
plot(u.sample, type="l")

