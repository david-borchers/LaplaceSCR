library(secr)
library(mvtnorm)

source("SCRutilities.r") # David's utility functions for SCR
library(fields) # needed for utility functions
library(corrplot) # for (scaled) second derivative matrix and correlation matrix plots

#--------------------------- Functions -------------------------------------------

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

# Calculates log(\prod_t\prod_k g_{itk}(s)) for a single individual for all s in mesh
# wi is ith row of binary proximity capture history (ith individual's capture history)
binprox.log.gi.si = function(wi,log.gk,log.gk1) {
  delta=rep(0,dim(log.gk)[1])
  delta[wi[wi>0]]=1
  log.gi.si <- delta %*% log.gk  + (1-delta) %*% log.gk1 # log(\prod_k g_itk(s))
  return(log.gi.si)
}


# Calculates log(\prod_t\prod_k g_{itk}(s)) for all individuals for all s in mesh
# capthist is n by nt binary proximity capture history matrix
# log.gk and log.gk1 are K by M matrices
# returns an n by M matrix
binprox.log.gi.s = function(capthist,log.gk,log.gk1) {
  log.gi.s=apply(capthist, 1, binprox.log.gi.si,log.gk=log.gk,log.gk1=log.gk1)
  return(t(log.gi.s))
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

# Probabilities for each detector and each mask point
# gk, log.gk and log.gk1 are all K x M matrices
gk <- g0 * exp(-dist^2 / 2 / sigma^2)
log.gk <- log(gk) # for use below
log.gk1 <- log(1-gk) # for use below

# probability of being caught at least once if at mesh vertex s of M
p..s <- 1 - apply(1-gk, 2, prod) ^ nt  # vector of length M

# calculate likelihoood components:

# Lambda.tilde
Lambda.tilde = sum(alpha*lambda*p..s);Lambda.tilde

# calculate G_i(s)
log.gi.s=binprox.log.gi.s(capthist,log.gk,log.gk1);dim(log.gi.s)
Gi.s=exp(log.gi.s)
Pi.s=t(alpha*lambda*t(Gi.s)) # transposing to multiply by row, not colum
Pi=apply(Pi.s,1,sum);Pi

# log-likelihood and likelihood:
loglik = -Lambda.tilde + sum(log(Pi)); loglik
lik=exp(loglik); lik

# First derivatives w.r.t. lambda
d.dlambda=(1/Pi) %*% Pi.s - alpha*lambda*p..s
#d.dlambda2=(1/Pi) %*% (alpha*lambda*exp(log.gi.s)) - alpha*lambda*p..s # (this is not correct)
#d.dlambda=alpha*lambda*((1/Pi) %*% exp(log.gi.s)) - alpha*lambda*p..s # (another way of calculating)
# Plot the first derivatives:
covariates(GMRFmesh)$d.dlambda=as.vector(d.dlambda)
plotcovariate(GMRFmesh,covariate="d.dlambda",main="First derivatives")
plot(dets,add=TRUE)
#plot(capthist,type="n.per.detector",add=TRUE)
plot(capthist,tracks=TRUE,rad=2,cappar=list(cex=0.75),add=TRUE)


# Second derivatives w.r.t. lambda
d2.dlambda2.k = alpha*lambda*((1/Pi) %*% exp(log.gi.s))^2
d2.dlambda2.j = -alpha*lambda
d2.dlambda2 = as.matrix(d2.dlambda2.j,nrow=1,drop=FALSE) %*% as.matrix(d2.dlambda2.k,ncol=1,drop=FALSE)
dim(d2.dlambda2)
diag(d2.dlambda2)=diag(d2.dlambda2) + d.dlambda

# Look at structure of 2nd derivative matrix 
# Warning: if M not very small, this is not very useful - difficult to make out
quartz()
div=max(abs(range(d2.dlambda))) # scale so in [-1, 1], so can use corrplot to display
corrplot(d2.dlambda/div,method="ellipse", main="scaled 2nd derivatives", cex.main=0.75)

# Calculate first derivatives of GMRF
f.e = dmvnorm(e, mean=rep(0,M), sigma=Sigma.e)
tau.e=solve(Sigma.e)
#quartz()
#div=max(abs(range(tau.e))) # scale so in [-1, 1], so can use corrplot to display
#corrplot(tau.e/div,method="ellipse", main="scaled precision", cex.main=0.75)
d.de = -f.e * tau.e %*% e
# plot for info:
covariates(GMRFmesh)$d.de=d.de
plotcovariate(GMRFmesh,covariate="d.de",main="GMRF first derivatives",contour=FALSE)
plot(dets,add=TRUE)

# Calculate second derivatives for GMRF
d2.de2 = f.e * tau.e %*% e %*% t(e) %*% tau.e - tau.e
# plot for info:
quartz()
div=max(abs(range(d2.de2))) # scale so in [-1, 1], so can use corrplot to display
corrplot(d2.de2/div,method="ellipse", main="scaled GMRF 2nd derivatives", cex.main=0.75)
hist(d2.de2)
det(d2.de2)
# Hmm ... machine preision problems? Else an error?

# Hessian:
d.dlambda.d.de = t(d.dlambda) %*% t(d.de)
H = d2.dlambda2*f.e + d.dlambda.d.de + t(d.dlambda.d.de) + lik*d2.de2
det(H)
# ... well that did not work! Hopefully machine precision issue?
