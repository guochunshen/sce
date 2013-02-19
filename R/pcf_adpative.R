#' estimate pair correlation function by adaptive method
#' 
#' \code{pcf_adaptive} is only a warp function of \link{mygest} function in InhomCluster package to
#' fit the formate of fv object defined in spatstat package. See details of the parameters in \link{mygest}.
#'  
#'  @param data a \code{spatstat ppp} point patttern object
#'  @param lambda vector containing intensity for each point in data
#'  @param maxt maximal distance t for which g(t) is estimated
#'  @param bw band with for kernel estimation
#'  @param adaptive if true, Snethlage (2001)'s adaptive estimate is used
#'  @param kerneltype if 1, the uniform kernel is used, 
#'  
#'  @details This function computes a kernel estimate of the g-function using the translation invariant edge correction (Ohser, 1983). The procedure only works for point patterns observed within a rectangular observation window.
#'  
#'  @return an \code{fv} object
#'  
#'  
pcf_adaptive<-function(data.ppp,lambda,maxt,bw,adaptive=T,kerneltype=1){
  g=mygest(data.ppp,lambda,maxt,bw,adaptive,kerneltype)

  g=data.frame(r=g$t,theo=rep(1,length(g$t)),trans=g$g)
  class(g)=c("fv","data.frame")
  attr(g,"argu")="r"
  attr(g,"valu")="trans"
  attr(g,"fmla")=". ~ r"
  attr(g,"alim")=range(g$r)
  attr(g,"units")="unit / units"
  attr(g,"labl")=c("r","{%s^{pois}}(r)","hat(%s^{trans})(r)")
  
  return(g)
}

mygest<-function (data, lambda, maxt, bw, adaptive = 1, kerneltype = 1) 
{
  if (data$window$type != "rectangle") {
    print("stop only rectangular observation windows allowed")
    stop
  }
  data$x = data$x - data$window$xrange[1]
  data$y = data$y - data$window$yrange[1]
  data$window$xrange = c(0, data$window$xrange[2] - data$window$xrange[1])
  data$window$yrange = c(0, data$window$yrange[2] - data$window$yrange[1])
  x = as.double(data$x)
  y = as.double(data$y)
  sideEW <- as.double(data$window$xrange[2])
  sideNS <- as.double(data$window$yrange[2])
  par = as.double(lambda)
  antpktdat = as.integer(length(x))
  maxt = as.double(maxt)
  bw = as.double(bw)
  gestimate = as.double(rep(0, 100))
  parametric = as.integer(F)
  correctype = as.integer(3)
  adaptive = as.integer(adaptive)
  kerneltype = as.integer(kerneltype)
  gestimate <- .C("gtrans", antpktdat, x, y, par, maxt, gestimate, 
                  bw, correctype, sideEW, sideNS, parametric, adaptive, 
                  kerneltype)[[6]]
  return(list(t = c(1:100) * maxt/100, g = gestimate))
}