## ----setup, include=FALSE------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- eval=FALSE---------------------------------------------------------
## install.packages("INLA", repos=c(getOption("repos"),
##         INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)


## ---- eval=FALSE---------------------------------------------------------
## library(INLA)
## inla.update(testing=T)


## ------------------------------------------------------------------------
library(INLA)
inla(formula=y~1,data=list(y=1:9))

