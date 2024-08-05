## ----setup, include=FALSE-----------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- eval=FALSE--------------------------------------------------
## ## Chunk not ran
## ## Example is real, but folder names were changed
## if (Sys.info()["user"]=="haakon") {
##   folder.data = "BarnacleDamage/"
##   folder.out = "../../../Barnacle/out-10-09/"
## } else if (Sys.info()["user"]=="tmjo") {
##   folder.data = "C:/Users/tmjo/Desktop/Haakon"
##   folder.out = "C:/Users/tmjo/Desktop/Haakon"
## } else {
##   stop("Add your folders as above.")
## }
## dir.create(folder.out, recursive = T, showWarnings = F)


## ---- eval=FALSE--------------------------------------------------
## ## Chunk not run
## ## Loading data
## filename = "data1"
## load(paste0(folder.data, filename))
## ## Writing a file
## filename = "figt1"
## png(file=paste0(folder.out, filename))
## # ...


## -----------------------------------------------------------------
## Velocity of the ship
v.ship = 10


## -----------------------------------------------------------------
n.row = 14
n.column = 19
for (i.row in 1:n.row) {
  for (i.column in 1:n.column) {
    ## code...
  }
}


## ---- eval=FALSE--------------------------------------------------
## ## Chunk not run
## ## Not:
## df
## df.scaled
## df2
## df2.scaled
## 
## ## Yes
## df1.raw
## df1.scaled
## df2.raw
## df2.scaled


## -----------------------------------------------------------------
y ~ elevation + aspect + dryness + vegetation + class + slope + ...


## ---- eval=FALSE--------------------------------------------------
## ## Chunk not run
## covar1.names = c("elevation", "aspect", "dryness", "vegetation", "class", "slope")
## Xcov1 = as.matrix(df1[, c("covar1.names")])
## y ~ Xcov1

