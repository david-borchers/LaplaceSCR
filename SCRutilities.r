
#' @title Plots image (and optionally contours) of mask covariate value
#' @param mask is an object of class `mask'
#' @param covariate is a character variable with the name of one of the covariates in mask (the one to plot)
#' @param contour is a logical, TRUE if want contour plots on image
#' @param ... other arguments to be passed to \code{prep4image}
plotcovariate=function(mask,covariate,contour=TRUE,key=TRUE, ...) {
  cnum=which(names(covariates(mask))==covariate)
  if(is.null(cnum)) stop(paste("No covariate(s) called",covariate))
  if(length(cnum)>1) warning("Can only plot one covariate at a time. First covariate being plotted.")
  dat=data.frame(x=mask$x,y=mask$y,z=covariates(mask)[[cnum]])
  prep4image(dat,contour=contour,key=key,...)
}

#' @title Prepares data frame for plotting with image/contour/persp.
#'   
#' @description From an input data frame with columns x, y and z, this function 
#'   creates a list with elements x, y and z in a format suitable for passing to
#'   functions \code{\link{image}}, \code{\link{contour}} or 
#'   \code{\link{persp}}. The coordinates in \code{data} are assumed to come
#'   from a 2D grid of points.
#'   
#' @param data a data frame with columns x and y being Cartesian coordinates, 
#'   and z being the values of some variable at each coordinate.
#' @param plot if \code{TRUE} then an image plot will be drawn using 
#'   \code{\link{image.plot}}
#' @param contour if \code{TRUE} then contours will be added (only used when 
#'   \code{plot=TRUE})
#' @param key logical for whether or not to include key when \code{plot = TRUE} (\code{\link{image.plot}} is used when \code{key = TRUE}, \code{\link{image}} is used when \code{key = FALSE})
#' @param ... other arguments to pass to \code{\link{image}} or \code{\link{image.plot}} (only used 
#'   when \code{plot=TRUE})
#'   
#' @details Sorts z on values of x first, then y, then creates a matrix of 
#'   z-values from this. Returns a list with elements x (unique values of x, in 
#'   increasing order), y (unique values of y, in increasing order) and z 
#'   (matrix of z-values in appropriate order for image/contour/persp). 
#'   
#'   If the original z is a factor variabele, the z returned is a matrix of integers 
#'   between 1 and length(levels(z)) and the output list has an attributes called 
#'   ``facnames'' that is a character vector containing the levels as factor 
#'   variables, with z=1 corresponding to the first name, z=2 to the second, etc.
#' @export
#' @importFrom fields image.plot
prep4image = function(data, plot = TRUE, contour = TRUE, key = TRUE, ...){
  
  # convert factor data$z to integer:
  zfactor=FALSE
  if(is.factor(data$z)) {
    zfactor=TRUE
    fac=data$z
    facnames=levels(fac)
    nlevels=length(facnames)
    data$z=rep(0,length(fac))
    got=rep(FALSE,nlevels)
    for(i in 1:nlevels){
      j=which(fac==facnames[i])
      if(length(j)>0) got[i]=TRUE
      data$z[j]=i
    }
    facnames=facnames[got] # remove factor names not in mask
  }
  data = as.matrix(data)
  
  x = sort(unique(data[,"x"]))
  y = sort(unique(data[,"y"]))
  
  z = matrix(NA, nr = length(x), nc = length(y))
  
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      m = which(data[,"x"] == x[i] & data[,"y"] == y[j]) ; m
      z[i,j] = if(length(m) == 0) NA else data[,"z"][m]
    }
  }
  
  if(plot){
    if(key){
      image.plot(x, y, z, ...)
    }else{
      image(x, y, z, ...)
    }
    if(contour) contour(x, y, z, add = TRUE)
  }
  
  outlist=list(x = x, y = y, z = z)
  if(zfactor) attributes(outlist)$facnames=facnames
  
  invisible(outlist)
  
}


#' @title add new covariate consisting of distances to points that have covariate==distance.to
#' @param mask is object of class "mask".
#' @param covariate is name of one of the elements of covariates(mask)
#' @param distance.to is the target value of covariate; distances to mask points with this value are calculated.
#' @param dname is the name you want to give to the new covariate
add.dist.to=function(mask,covariate,distance.to,dname=NULL,overwrite=FALSE){
  covs=names(covariates(mask))
  if(is.null(covs)) stop("No covariates in mask. Can't add distance to anything")
  ncov=length(covs)
  if(is.null(dname)) dname=paste("dist2",covariate,"=",distance.to,sep="")
  if(is.element(covariate,covs)) {
    if(dname %in% names(covariates(mask))) {
      if(overwrite) {
        warning(paste("Covariate called",dname,"already existed. It has been overwritten."))
        covi=which(names(covariates(mask))==dname)
      } else stop(paste("Covariate called",covariate,"already exists. Use overwrite=TRUE if you want to overwrite it."))
    } else {
      covi=ncov+1
    }
    cov=covariates(mask)[which(covs==covariate)]
    if(is.null(cov)) stop(paste("No covariate called",covariate,"in mask."))
    targets=which(cov==distance.to)
    distances=rdist(mask,mask[targets,])
    d2=apply(distances,1,min)
    covariates(mask)[[covi]]=d2
    names(covariates(mask))[[covi]]=dname
    return(mask)
  } else stop(paste("No covariate called",covariate,"in mask."))
}


