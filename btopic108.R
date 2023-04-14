## ----setup, include=FALSE-----------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE, message=FALSE--------------------------------
library(INLA); library(sp); library(fields)
library(geoR)
library(viridisLite)
# - for better colours
rm(list=ls())
options(width=70, digits=2)


## -----------------------------------------------------------------
data('ca20')
class(ca20)
summary(ca20)


## -----------------------------------------------------------------
df = data.frame(y = ca20$data, locx = ca20[[1]][ , 1], locy = ca20[[1]][ , 2], ca20[[3]])
spatial.scaling = 100
df$locx = (df$locx - min(df$locx))/spatial.scaling
df$locy = (df$locy - min(df$locy))/spatial.scaling
df$altitude = df$altitude - mean(df$altitude)
df$y = df$y-50


## -----------------------------------------------------------------
summary(df)

head(df)

cor(cbind(df[, 1:4], as.numeric(df[ , 5])))


## -----------------------------------------------------------------
plot(df$altitude, df$y)
abline(lm(df$y~df$altitude), col="red")


## -----------------------------------------------------------------
quilt.plot(x=df$locx,y=df$locy,z=df$y,nx=40,ny=40, col = plasma(101),
           main = "Data")


## -----------------------------------------------------------------
max.edge = 0.5
mesh <- inla.mesh.2d(
  loc=df[ , c('locx', 'locy')],
  offset = c(0.5, 1.5),
  max.edge=c(max.edge, max.edge*3),
  # discretization accuracy
  cutoff=max.edge/5)
# cutoff removes locations that are too close, good to have >0


## -----------------------------------------------------------------
plot(mesh, asp=1)
points(df[ , c('locx', 'locy')], col='red')
axis(1); axis(2)


## -----------------------------------------------------------------
A = inla.spde.make.A(mesh=mesh, loc=data.matrix(df[ , c('locx', 'locy')]))
dim(A)
A[1:2, 100:200]
# - to see what A represents


## -----------------------------------------------------------------
Xcov = data.frame(intercept=1, altitude=df$altitude)
# - could add: area1 = (df$area==1)*1, area2 = (df$area==2)*1
# - - expands the factor covariates
# - ensure that all entries are numeric!
Xcov = as.matrix(Xcov)
colnames(Xcov)

stack <- inla.stack(tag='est',
                    # - Name (nametag) of the stack
                    # - Here: est for estimating
                    data=list(y=df$y),
                    effects=list(
                      # - The Model Components
                      s=1:mesh$n,
                      # - The "s" is means "spatial"
                     Xcov=Xcov),
                      # - The second is all fixed effects
                    A = list(A, 1)
                    # - First projector matrix is for 's'
                    # - second is for 'fixed effects'
                    )


## -----------------------------------------------------------------
prior.median.sd = 1; prior.median.range = 7
# - diff(range(mesh$loc[, 1]))/2
# - sd(df$y)/10
# - thisk about these, and experiment!
spde = inla.spde2.pcmatern(mesh, prior.range = c(prior.median.range, .5), prior.sigma = c(prior.median.sd, .5), constr = T)


## -----------------------------------------------------------------
formula = y ~ -1 + Xcov + f(s, model=spde)
# - Remove standard intercept
# - Fixed effects + random effects


## -----------------------------------------------------------------
prior.median.gaus.sd = 5.5
# - Think about this value
# - Remember sd(df$y)
family = 'gaussian'
control.family = list(hyper = list(prec = list(
  prior = "pc.prec", fixed = FALSE, param = c(prior.median.gaus.sd,0.5))))


## -----------------------------------------------------------------
res <- inla(formula, data=inla.stack.data(stack),
            control.predictor=list(A = inla.stack.A(stack), compute=T),
            # compute=T to get posterior for fitted values
            family = family,
            control.family = control.family,
            #control.compute = list(config=T, dic=T, cpo=T, waic=T), 
            # - Model comparisons
            control.inla = list(int.strategy='eb'),
            # - faster computation
            #control.inla = list(int.strategy='grid'),
            # - More accurate integration over hyper-parameters
            verbose=F)


## -----------------------------------------------------------------
summary(res)


## -----------------------------------------------------------------
for (i in 1:length(res$marginals.fixed)) {
  tmp = inla.tmarginal(function(x) x, res$marginals.fixed[[i]]) 
  plot(tmp, type = "l", xlab = paste("Fixed effect marginal", i, ":", res$names.fixed[i]), ylab = "Density")
}


## -----------------------------------------------------------------
tmp = inla.tmarginal(function(x) exp(-x), res$internal.marginals.hyperpar[[2]]) 
plot(tmp, type = "l", xlab = "inverse range", ylab = "Density")
xvals = seq(0, 10, length.out=1000)
lambda = -log(.5)/(1/prior.median.range); lines(xvals, 6*exp(-lambda*xvals), lty='dashed')


## -----------------------------------------------------------------
tmp = inla.tmarginal(function(x) exp(x), res$internal.marginals.hyperpar[[3]]) 
plot(tmp, type = "l", xlab = expression(sigma[u]), ylab = "Density")
xvals = seq(1, 20, length.out=1000)
lambda = -log(.5)/prior.median.sd; lines(xvals, 20*exp(-lambda*xvals), lty='dashed')


## -----------------------------------------------------------------
tmp = inla.tmarginal(function(x) exp(-0.5*x), res$internal.marginals.hyperpar[[1]]) 
plot(tmp, ty = "l", xlab = expression(sigma[iid]), yla = "Density")
xvals = seq(0, 10, length.out=1000)
lambda = -log(.5)/prior.median.gaus.sd; lines(xvals, .5*exp(-lambda*xvals), lty='dashed')


## -----------------------------------------------------------------
local.plot.field = function(field, mesh, xlim=c(0,11), ylim=c(0,9), ...){
  stopifnot(length(field) == mesh$n)
  # - error when using the wrong mesh
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 plotting grid 
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, col = plasma(101), ...)  
}


## -----------------------------------------------------------------
local.plot.field(res$summary.random[['s']][['mean']], mesh)
lines(5+c(-0.5, 0.5)*(res$summary.hyperpar[2, '0.5quant']), c(1,1)*5, lwd=3)
# - add on the estimated range
axis(1); axis(2)
# - the transformed axes
# - could have used the original scale


## -----------------------------------------------------------------
local.plot.field(res$summary.random$s$sd, mesh)
axis(1); axis(2)


## -----------------------------------------------------------------
quilt.plot(x=df$locx,y=df$locy,z=res$summary.fitted.values$mean[1:nrow(df)],nx=40,ny=40, col = plasma(101), main="Fitted values", 
           zlim = range(df$y))

