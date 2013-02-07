#' estimate pair correlation function by adaptive method
#' 
#' \code{pcf_adaptive} is only a warp function of \link{mygest} function in InhomCluster package to
#' fit the formate of fv object defined in spatstat package. See details of the parameters in \link{mygest}.
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