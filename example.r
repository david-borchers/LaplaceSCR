library(secr)
library(mvtnorm)

source("SCRutilities.r") # David's utility functions for SCR
library(fields) # needed for utility functions

# set detector type
detype="multi"
detype="count"

# create the detectors
spacing=20
dets = make.grid(nx=7, ny=7, spacing=spacing,detector=detype)
plot(dets,border=0)
# set normal detection function parameters
g0=0.5; sigma=spacing
# create the mesh
mesh0 = make.mask(dets, buffer=4*sigma, nx=20, ny=20, type="trapbuffer")

# make GMRF with speficied variance and covariance range to mesh, such that E(N) about 100
gmrf=make.gmrf(mesh0,sigma.xi=1,covrange=25,E.N=200,seed=1246)
#gmrf=make.gmrf(mesh0,sigma.xi=1,covrange=300,E.N=100)
# Look at intensity surface
mesh=mesh0
covariates(mesh)$lambda = gmrf$lambda
plotcovariate(mesh,covariate="lambda",contour=FALSE)
plot(dets,add=TRUE) # ... with detectors overlaid
gmrf$E.N
summary(gmrf$lp)
xi = gmrf$xi
Sigma.xi = gmrf$Sigma.xi
# Look at covariance matrix
M=dim(mesh)[1]
image.plot(1:M,1:M, gmrf$Sigma.xi, main="GMRF Variance-covariance")

# simulate a population from this intensity surface
pop=sim.popn(D="lambda", core=mesh, model2D="IHP", seed=12345)
# plot mesh with individuals' locations and detectors overlaid
plot(mesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets,add=TRUE)
dim(pop)[1] # check simulated population size

# Generate capture histories
# first set number of occasions
capthist=sim.capthist(dets,popn=pop, detectpar=list(g0=g0,sigma=sigma), noccasions=nt, nsessions=1,seed=12345)
summary(capthist)
n=dim(capthist)[1]
plot(mesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets,add=TRUE)
plot(capthist, border=sigma, tracks=TRUE, gridlines=FALSE,rad=3,add=TRUE)

# make list of SCR stuff to pass
scr.list = make.scr.list(capthist, mesh, g0, sigma, base.lp=gmrf$lp)

# Calculate the log-likelihood:
# =============================
ll = scr.loglik(xi, Sigma.xi, scr.list)
ll$loglik;ll$scr;ll$gmrf

# For interest, look at where joint scr and gmrf likelihood is maximised
# searching over a range of multiples of xi:
nxi=10
fac=seq(0,2,length=nxi)
ll.tot=ll.scr=ll.gmrf=rep(NA,nxi)
for(i in 1:nxi) {
  ll = scr.loglik(xi*fac[i], Sigma.xi, scr.list)
  ll.tot[i] = ll$loglik
  ll.scr[i] = ll$scr
  ll.gmrf[i] = ll$gmrf
}
plot(fac,ll.tot,type="l",main="joint likelihood")
plot(fac,ll.scr,type="l",main="scr component of joint likelihood")
plot(fac,ll.gmrf,type="l",main="gmrf component of joint likelihood") # (must always have max at xi=0)



# Calculate the Hessian and other derivatives:
# ===========================================
derivs=scr.derivs(xi, Sigma.xi, scr.list)

# Plot the first derivatives of l:
covariates(mesh)$first.l=as.vector(derivs$first.l)
plotcovariate(mesh,covariate="first.l",main="First derivatives",contour=FALSE)
plot(dets,add=TRUE) # add detectors
#plot(capthist,tracks=TRUE,rad=2,cappar=list(cex=0.75),add=TRUE) # ... and detections

# Look at diagonal of 2nd derivatives of l in space:
covariates(mesh)$second.l=diag(derivs$second.l)
plotcovariate(mesh,covariate="second.l",main="Diagonal 2nd derivatives of l",contour=FALSE)
plot(dets,add=TRUE)

# Look at structure of 2nd derivative matrix of log.l
image.plot(1:M,1:M,derivs$second.l, main="2nd derivatives of l")
det(derivs$second.l)
# Zero determinant not so good! ... machine precision problems? Else an error?

# Plot the first derivatives of log GMRF:
covariates(mesh)$first.loggmrf=as.vector(derivs$first.loggmrf)
plotcovariate(mesh,covariate="first.loggmrf",main="First derivatives of log GMRF",contour=FALSE)
plot(dets,add=TRUE)
#plot(capthist,tracks=TRUE,rad=2,cappar=list(cex=0.75),add=TRUE)

# Look at 2nd derivative matrix of log GMRF
image.plot(1:M,1:M,derivs$second.loggmrf, main="2nd derivatives of log GMRF")
det(solve(derivs$second.loggmrf))
# Zero determinant not so good

# Look at diagonal of 2nd derivatives of log GMRF in space:
covariates(mesh)$second.loggmrf=diag(derivs$second.loggmrf)
plotcovariate(mesh,covariate="second.loggmrf",main="GMRF diagnonal 2nd derivatives",contour=FALSE)
plot(dets,add=TRUE)

# Look at Hessian
image.plot(1:M,1:M,derivs$H, main="Hessian")
det(solve(derivs$H))
# Zero determinant still not so good.

