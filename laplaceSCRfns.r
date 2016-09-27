#' @title First and second derivatives of SCR log-likelihood w.r.t. GMRF intensities
#'
#' @description  Calculates the first and second derivatives of SCR log-likelihood hwith respect to
#' the GMRF intensity at every point in the mesh.
#'
#' @param alpha The integration weight attachced to each mesh point: a vector of length M
#' @param lambda The Poisson intensity at each mesh point: a vector of length M
#' @param gtk A KxM matrix with each row corresponding to a detector and each column to a mesh point 
#' (in the same order as the rows in mesh)
#' @param xi values of GMRF at each mesh point: Mx2 matrix with x- and y- coordinates of mesh points in
#' first and second columns
#' @param Sigma.xi MxM variance-covariance matrix of xi.
#' @param capthist A capture history object in the form required by package \code{secr}. Currently only 
#' ONE SESSION capture histories are supported.
#' 
#' @return A list with each element being an MxM matrix of derivtives: 
#' \itemize{
#'  \item{"H"}{Hessian of marginal likelihood}
#'  \item{"first"}{First derivatives of marginal likelihood}
#'  \item{"first.l"}{first derivative of likelihood given GMRF}
#'  \item{"second.l"}{second derivative of likelihood given GMRF}
#'  \item{"first.loggmrf"}{first derivative of log GMRF}
#'  \item{"second.loggmrf"}{second derivative of log GMRF}
#'  }
#' 
#' @references David and Finn's document ``SCR Likelihood with a spatial log Gaussian Cox process and 
#' GMRF''.
#'
#' @examples
#' 
#' @export
scr.derivs=function(xi, Sigma.xi, scr.list) {
  
  lambda = exp(scr.list$base.lp + xi)
  alpha = scr.list$alpha

  # alpha_j * lambda_j * p..(s_j): vector of length M
  alp..s = alpha * lambda * scr.list$p..s

  # alpha_j * lambda_j * P_i(s_j): nxM matrix
  alPi.s = t(alpha * lambda * t(scr.list$Pi.s)) # transposing to multiply by row, not colum
  # sum_j (alpha_j * lambda_j * P_i(s_j)): vector of length n
  alPi = apply(alPi.s,1,sum)
  # [alpha_j * lambda_j * P_i(s_j)] / [sum_j (alpha_j * lambda_j * P_i(s_j))]: vector of length M
  alPisum = (1/alPi) %*% alPi.s
  
  # First derivatives w.r.t. lambda
  dl.dlambda = alPisum - alp..s

  # Second derivatives w.r.t. lambda
  d2l.dlambda2 = -t(alPisum) %*% alPisum
  diag(d2l.dlambda2) = diag(d2l.dlambda2) + dl.dlambda
  
  # Calculate first derivatives of log GMRF wrt xi
  tau.xi = solve(Sigma.xi)
  dlogf.dxi = tau.xi %*% xi
  
  # Calculate second derivatives for log GMRF wrt log(xi)
  d2logf.dxi2 = tau.xi
  
  # Hessian:
  H = d2l.dlambda2 + d2logf.dxi2
  
  return(list(H=H, first=dl.dlambda+t(dlogf.dxi), first.l=dl.dlambda, second.l=d2l.dlambda2, 
              first.loggmrf=dlogf.dxi, second.loggmrf=d2logf.dxi2))
}


#' @title Joint SCR log-likelihood with GMRF
#'
#' @description  Calculates the joint SCR and GMRF log-likelihood.
#'
#' @param alpha The integration weight attachced to each mesh point: a vector of length M
#' @param lambda The Poisson intensity at each mesh point: a vector of length M
#' @param gtk A KxM matrix with each row corresponding to a detector and each column to a mesh point 
#' (in the same order as the rows in mesh)
#' @param xi values of GMRF at each mesh point: Mx2 matrix with x- and y- coordinates of mesh points in
#' first and second columns
#' @param Sigma.xi MxM variance-covariance matrix of xi.
#' @param capthist A capture history object in the form required by package \code{secr}. Currently only 
#' ONE SESSION capture histories are supported.
#' 
#' @return Value of log-likelihood
#' 
#' @references David and Finn's document ``SCR Likelihood with a spatial log Gaussian Cox process and 
#' GMRF''.
#'
#' @examples
#' 
#' @export
scr.loglik=function(xi, Sigma.xi, scr.list) {

  lambda = exp(scr.list$base.lp + xi)
  alp..s = scr.list$alpha * lambda * scr.list$p..s
  Lambda.tilde = sum(alp..s)
  
  # calculate alpha_j * lambda_j * P_i(s_j)
  alPi.s=t(scr.list$alpha*lambda*t(scr.list$Pi.s)) # transposing to multiply by row, not colum
  alPi=apply(alPi.s,1,sum)
  
  # log-likelihood given GMRF
  l = sum(log(alPi)) - Lambda.tilde
  
  # GMRF component
  log.f.xi = log(dmvnorm(xi, mean=rep(0,M), sigma=Sigma.xi))
  
  loglik = l + log.f.xi
  return(list(loglik=loglik, scr=l, gmrf=log.f.xi))
}


make.scr.list=function(capthist, mesh, g0, sigma, base.lp=0) {
  # Calculate some things needed in likelihood and Hessian that do not depend on lambda
  # -----------------------------------------------------------------------------------
  M=dim(mesh)[1] # number of mesh points
  alpha=rep(1/M,M) # integration weights
  
  # calculate distances from each detector to each mesh point
  dist=distances(traps(capthist),mesh)

  # Proximity detector with count data
  if (detector(traps(capthist))=="count") {
    # Probabilities for each detector and each mesh point (assumed same at all times)
    er <- g0 * exp(-dist^2 / 2 / sigma^2) # Poisson encounter rate at all M distances from each detector
    
    # probability of being caught at least once if at mesh vertex s of M
    p..s <- 1 - exp(- apply(er * dim(capthist)[2],2,sum))
    
    # calculate P_i(s): nxM matrix
    comb.capthist=apply(capthist,c(1,3),sum)
    log.Pi.s=count.log.Pi.s(comb.capthist,er)
    Pi.s=exp(log.Pi.s)
  }
  # Multi-catch trap
  if (detector(traps(capthist))=="multi") {
    # Probabilities for each detector and each mesh point (assumed same at all times)
    gtk <- g0 * exp(-dist^2 / 2 / sigma^2)
    log.gtk <- log(gtk) # for use below
    log.gtk1 <- log(1-gtk) # for use below
    
    # probability of being caught at least once if at mesh vertex s of M
    p..s <- 1 - apply(1-gtk, 2, prod) ^ dim(capthist)[2] 
    
    # calculate P_i(s): nxM matrix
    log.Pi.s=multi.log.Pi.s(capthist,log.gtk,log.gtk1)
    Pi.s=exp(log.Pi.s)
  }
  
  # put all this crap in a list to pass
  return(list(alpha=alpha, p..s=p..s, Pi.s=Pi.s, base.lp=base.lp))
}





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

# MULTI-CATCH detectors (detector(traps(capthist))=="multi")
# ---------------------
# Calculates log(P_i(s)) for a single individual for all s in mesh
# wi is ith row of multi-catch trap capture history (ith individual's capture history)
# of length nt, and with detector number in non-zero elements
multi.log.Pi.si = function(wi,log.gtk,log.gtk1) {
  delta=rep(0,dim(log.gtk)[1])
  delta[wi[wi>0]]=1
  log.Pi.si <- delta %*% log.gtk  + (1-delta) %*% log.gtk1 # log(\prod_k g_itk(s))
  return(log.Pi.si)
}

# Calculates log(P_i(s)) for all individuals for all s in mesh
# capthist is n by nt multi-catch trap capture history matrix
# log.gtk and log.gtk1 are K by M matrices
# returns an n by M matrix
multi.log.Pi.s = function(capthist,log.gtk,log.gtk1) {
  log.Pi.s=apply(capthist, 1, multi.log.Pi.si,log.gtk=log.gtk,log.gtk1=log.gtk1)
  return(t(log.Pi.s))
}


# COUNT detectors (detector(traps(capthist))=="count")
# ---------------------
# Calculates log(P_i(s)) for a single individual for all s in mesh
# wi is ith row of count capture history (ith individual's capture history)
# of length nt, and with detector number in non-zero elements
count.log.Pi.si = function(wi,er) {
  one=rep(1,length(wi))
  log.Pi.si = wi %*% log(er)  - one %*% er - sum(lfactorial(wi)) # log(\prod_k Poisson(n_{ks}))
  return(log.Pi.si)
}

# Calculates log(P_i(s)) for all individuals for all s in mesh
# comb.capthist is n by K matrix of counts ("count" capthist collapsed over occasions)
# log.gtk and log.gtk1 are K by M matrices
# returns an n by M matrix
count.log.Pi.s = function(comb.capthist,er) {
  log.Pi.s=apply(comb.capthist, 1, count.log.Pi.si,er=er)
  return(t(log.Pi.s))
}

#' @title Generate GMRF and intensity on a mesh
#'
#' @description  Generates a GMRF and intensity for every point on a mesh.
#'
#' @param mesh A 'mask', as required by package \code{secr}
#' @param sigma.xi  Variance of GMRF at distance zero
#' @param covrange Range of Matern covariance function
#' @param seed Random number seed (for repeatability)
#' @param E.N target E(N) in region
#' @param center Logical. If TRUE, realised GMRF values are centred about zero, else not
#' 
#' @examples
#' library(secr)
#' # create the detectors
#' spacing=20
#' dets = make.grid(nx=5, ny=5, spacing=spacing)
#' plot(dets,border=0)
#' # set normal detection function parameters
#' g0=0.5; sigma=spacing
#' # create the mesh
#' mesh = make.mask(dets, buffer=4*sigma, nx=20, ny=20, type="trapbuffer")
#' 
#' # make GMRF with speficied variance and covariance range to mesh, such that E(N)=100
#' gmrf=make.gmrf(mesh,sigma.xi=1,covrange=300,E.N=100,seed=2016)
#' # Plot covariance matrix:
#' M=dim(mesh)[1]
#' image.plot(1:M,1:M, gmrf$Sigma.xi, main="GMRF Variance-covariance")
#' covariates(mesh)$lambda = gmrf$lambda
#' sum(covariates(mesh)$lambda*attributes(mesh)$area) # check this equals N
#' plotcovariate(mesh,covariate="lambda",contour=FALSE) # Look at it
#' plot(dets,add=TRUE) # ... with detectors overlaid
#' 
make.gmrf=function(mesh, sigma.xi, covrange, E.N=NULL, seed=NULL, center=TRUE) {
  # calculate distances between mesh points
  mdist=distances(mesh,mesh)
  
  # generate GMRF and add it and lambda to mesh covariates
  Corr.xi = apply(mdist,1,FUN="Matern",range=covrange) # Matern correlation matrix
  Sigma.xi=Corr.xi * sigma.xi^2 # GMRF Variance matrix
  meanxi=0
  mu.xi = rep(meanxi,dim(mesh)[1]) # GMRF mean
  # generate GMRF variates
  set.seed(seed=seed)
  xi = as.vector(rmvnorm(1,mean=mu.xi,sigma=Sigma.xi))
  if(center) {
    meanxi=mean(xi)
    xi=xi-meanxi
  }

  a = attributes(mesh)$area # area of each mesh cell
  if(is.null(E.N)) {
    if(!is.element("lambda",names(covariates(mesh)))) {
      lp=0 + meanxi # systematic part of linear predictor
    } else {
      lp=log(covariates(mesh)$lambda) + meanxi - (sigma.xi^2)/2 # systematic part of linear predictor
    }
  } else {
    lambda0 = E.N/(a*dim(mesh)[1])
    if(!is.element("lambda",names(covariates(mesh)))) {
      lp = log(lambda0) + meanxi - (sigma.xi^2)/2 # systematic part of linear predictor
    } else {
      lp = log(covariates(mesh)$lambda*lambda0/mean(covariates(mesh)$lambda)) + meanxi - (sigma.xi^2)/2 # systematic part of linear predictor
    }
  }

  lambda = exp(lp + xi); sum(lambda*a)
  
  return(list(lambda=lambda, xi=xi, lp=lp, Sigma.xi=Sigma.xi, approx.E.N=sum(lambda*a)))
}
