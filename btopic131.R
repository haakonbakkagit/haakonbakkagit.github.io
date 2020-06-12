## ----setup, include=FALSE-------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- eval=FALSE----------------------------------------------------
## ## Chunk not ran
## ## Example is real, but folder names were changed
## if (Sys.info()["user"]=="haakon") {
##   folder.data = "BarnacleDamage/"
##   folder.out = "../../../Barnacle/out-10-09/"
## } else if (Sys.info()["user"]=="tmjo") {
##   folder.data = "C:/Users/tmjo/Desktop/Haakon"
##   folder.out = "C:/Users/tmjo/Desktop/Haakon"
## } else {
##   stop("Add your folders as above:")
## }
## dir.create(folder.out)


## ---- eval=FALSE----------------------------------------------------
## ## Chunk not ran
## ## Loading data
## filename = "data1"
## load(paste0(folder.data, filename))
## ## Writing a file
## filename = "figt1"
## png(file=paste0(folder.out, filename))
## # ...

