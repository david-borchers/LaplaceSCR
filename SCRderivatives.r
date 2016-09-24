library(secr)
library(mvtnorm)

source("SCRutilities.r") # David's utility functions for SCR
library(fields) # needed for utility functions
library(corrplot) # for (scaled) second derivative matrix and correlation matrix plots

#--------------------------- Functions -------------------------------------------

# Calucates distances between two sets of coordinates
# X is a by 2 matrix, Y is b by 2 matrix, output is a by b matrix
distances <- function (X, Y) {
  ## X and Y are 2-column matrices of coordinates
  onerow <- function (xy) {
    d <- function(xy2) {
      sqrt(sum((xy2 - xy)^2))
    }
    apply(Y, 1, d)
  }
  t(apply(X, 1, onerow))
}

# Calculates log(P_i(s)) for a single individual for all s in mesh
# wi is ith row of binary proximity capture history (ith individual's capture history)
# of length nt, and with detector number in non-zero elements
binprox.log.Pi.si = function(wi,log.gtk,log.gtk1) {
  delta=rep(0,dim(log.gtk)[1])
  delta[wi[wi>0]]=1
  log.Pi.si <- delta %*% log.gtk  + (1-delta) %*% log.gtk1 # log(\prod_k g_itk(s))
  return(log.Pi.si)
}

# Calculates log(P_i(s)) for all individuals for all s in mesh
# capthist is n by nt binary proximity capture history matrix
# log.gtk and log.gtk1 are K by M matrices
# returns an n by M matrix
binprox.log.Pi.s = function(capthist,log.gtk,log.gtk1) {
  log.Pi.s=apply(capthist, 1, binprox.log.Pi.si,log.gtk=log.gtk,log.gtk1=log.gtk1)
  return(t(log.Pi.s))
}

#----------------------------------------------------------------------



# create the detectors
spacing=20
dets = make.grid(nx=5, ny=5, spacing=spacing)
plot(dets,border=0)

# set normal detection function parameters
g0=0.5; sigma=spacing

# create the mesh
mesh = make.mask(dets, buffer=4*sigma, nx=20, ny=20, type="trapbuffer")
plot(mesh, border=5,dots=FALSE, col="white", meshcol="gray")
plot(dets, add=TRUE)
M=dim(mesh)[1] # number of mesh points
alpha=rep(1/M,M) # integration weights


# Create a GMRF
mdist=distances(mesh,mesh)
dim(mdist)
sigma.e = 1 # Variance at distance zero
covrange=300 # range of Matern covariance function
#cov.e = Matern(mdist[1,], range=covrange) * sigma.e^2
#plot(cov.e)
Corr.e = apply(mdist,1,FUN="Matern",range=covrange) # Matern correlation matrix
#quartz()
#corrplot(Corr.e,method="ellipse", main="GMRF correlation", cex.main=0.75)
Sigma.e=Corr.e * sigma.e^2 # GMRF Variance matrix
mu.e = rep(0,M) # GMRF mean
# generate some GMRF variates
set.seed(seed=2016)
e = as.vector(rmvnorm(1,mean=mu.e,sigma=Sigma.e))
GMRFmesh=mesh
covariates(GMRFmesh)$e = e
lambda.lp.intercept = log(2)
covariates(GMRFmesh)$lambda = exp(lambda.lp.intercept + as.vector(e))
#plotcovariate(GMRFmesh,covariate="e")
plotcovariate(GMRFmesh,covariate="lambda")
#median(covariates(GMRFmesh)$lambda)
#mean(covariates(GMRFmesh)$lambda)
plot(dets,add=TRUE)

# Scale lambda to give E(N)=N
N=100
a=attributes(mesh)$area # area of each mesh cell
# add lambda to the mesh as a covariate
covariates(GMRFmesh)$lambda=N*covariates(GMRFmesh)$lambda/(sum(covariates(GMRFmesh)$lambda)*a)
lambda=covariates(GMRFmesh)$lambda # (so don't have to type it all out below)
sum(covariates(GMRFmesh)$lambda*a) # check this equals N
plotcovariate(GMRFmesh,covariate="lambda") # Look at it
plot(dets,add=TRUE) # ... with detectors overlaid

# simulate a population from this intensity surface
pop=sim.popn(D="lambda", core=GMRFmesh, model2D="IHP", seed=12345)
# plot mesh with individuals' locations and detectors overlaid
plot(GMRFmesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets,add=TRUE)
dim(pop)[1] # check simulated population size

# Generate capture histories
nt=5 # number of occasions
capthist=sim.capthist(dets,popn=pop, detectpar=list(g0=g0,sigma=sigma), noccasions=nt, nsessions=1,seed=12345)
summary(capthist)
n=dim(capthist)[1]
plot(GMRFmesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets,add=TRUE)
plot(capthist, border=sigma, tracks=TRUE, gridlines=FALSE,rad=3,,add=TRUE)

# calculate distances from each detector to each mesh point
dist=distances(dets,mesh)

# Probabilities for each detector and each mesh point (assumed same at all times)
# gk, log.gtk and log.gtk1 are all K x M matrices
gtk <- g0 * exp(-dist^2 / 2 / sigma^2)
log.gtk <- log(gtk) # for use below
log.gtk1 <- log(1-gtk) # for use below

# probability of being caught at least once if at mesh vertex s of M
# vector of length M
p..s <- 1 - apply(1-gtk, 2, prod) ^ nt  

# calculate likelihoood components:

# Lambda.tilde
Lambda.tilde = sum(alpha*lambda*p..s);Lambda.tilde

# calculate P_i(s)
log.Pi.s=binprox.log.Pi.s(capthist,log.gtk,log.gtk1);dim(log.Pi.s)
Pi.s=exp(log.Pi.s)
Pi.s=t(alpha*lambda*t(Pi.s)) # transposing to multiply by row, not colum
Pi=apply(Pi.s,1,sum);Pi

# First derivatives w.r.t. lambda
d.dlambda=(1/Pi) %*% Pi.s - alpha*lambda*p..s
#d.dlambda2=(1/Pi) %*% (alpha*lambda*Pi.s) - alpha*lambda*p..s # (this is not correct)
#d.dlambda=alpha*lambda*((1/Pi) %*% Pi.s) - alpha*lambda*p..s # (another way of calculating)
# Plot the first derivatives:
covariates(GMRFmesh)$d.dlambda=as.vector(d.dlambda)
plotcovariate(GMRFmesh,covariate="d.dlambda",main="First derivatives")
plot(dets,add=TRUE)
#plot(capthist,type="n.per.detector",add=TRUE)
plot(capthist,tracks=TRUE,rad=2,cappar=list(cex=0.75),add=TRUE)


# Second derivatives w.r.t. lambda
d2.dlambda2.k = alpha*lambda*((1/Pi) %*% Pi.s)
d2.dlambda2.j = -alpha*lambda*((1/Pi) %*% Pi.s)
#d2.dlambda2 = as.matrix(d2.dlambda2.j,nrow=1,drop=FALSE) %*% as.matrix(d2.dlambda2.k,ncol=1,drop=FALSE)
d2.dlambda2 = t(d2.dlambda2.j) %*% d2.dlambda2.k
dim(d2.dlambda2)
diag(d2.dlambda2)=diag(d2.dlambda2) + d.dlambda
det(d2.dlambda2)
# Zero determinant is not good!

# Try to look at structure of 2nd derivative matrix 
quartz()
#div=max(abs(range(d2.dlambda2))) # scale so in [-1, 1], so can use corrplot to display
#corrplot(d2.dlambda2/div,method="ellipse", main="scaled 2nd derivatives", cex.main=0.75)
image.plot(1:M,1:M,d2.dlambda2)


# Calculate first derivatives of log(xi)
#f.e = dmvnorm(e, mean=rep(0,M), sigma=Sigma.e)
tau.e=solve(Sigma.e)
#quartz()
#div=max(abs(range(tau.e))) # scale so in [-1, 1], so can use corrplot to display
#corrplot(tau.e/div,method="ellipse", main="scaled precision", cex.main=0.75)
#d.de = -f.e * tau.e %*% e
dlogf.de = tau.e %*% e
# plot for info:
covariates(GMRFmesh)$dlogf.de=dlogf.de
plotcovariate(GMRFmesh,covariate="dlogf.de",main="GMRF first derivatives",contour=FALSE)
plot(dets,add=TRUE)

# Calculate second derivatives for GMRF
#d2.de2 = f.e * tau.e %*% e %*% t(e) %*% tau.e - tau.e
d2logf.de2 = tau.e
# plot for info:
quartz()
#div=max(abs(range(d2logf.de2))) # scale so in [-1, 1], so can use corrplot to display
#corrplot(d2logf.de2/div,method="ellipse", main="scaled GMRF 2nd derivatives", cex.main=0.75)
image.plot(1:M,1:M,d2logf.de2)
hist(d2logf.de2)
det(d2logf.de2)
# Hmm ... machine preision problems? Else an error?

# Look at diagonal of 2nd derivative matrix in space:
covariates(GMRFmesh)$d2logf.de2=diag(d2logf.de2)
plotcovariate(GMRFmesh,covariate="d2logf.de2",main="GMRF diagnonal 2nd derivatives",contour=FALSE)
plot(dets,add=TRUE)

# Hessian:
H = d2.dlambda2 + d2logf.de2
image.plot(1:M,1:M,H)
det(H)
# ... well that determinant is not good either! Hopefully machine precision issue?
