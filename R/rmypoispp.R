#'
#' generate a poisson distribution point pattern, avoid 0 and very large number of point
#'
#' @param lambda an im object represent the intensity map of the expected point pattern
#' @param Nexpect the expect number of points in the observed window
#' 
#' @details
#' This is a warp function of \code{\link{rpoispp}} in spatstat package, and add two limitations into it. one is that avoid 
#' generating too much data points that commonly generated by a wrong lambda. thus it avoid the memory problem. The second 
#' improve is elimate the possibility of generating a point pattern with 0 point. this commonly happened in the case of very 
#' small lambda. although there is no problem with the zero point pattern, it will cause any other unexpected problem in other
#' functions that use the simulated data point.
#'
#'

rmypoispp<-function(lambda,Nexpect){
  lambda=check_lambda(lambda,Nexpect)
  data.ppp=rpoispp(lambda)
  while(data.ppp$n<1){
    data.ppp=rpoispp(lambda)
  }
  
  return(data.ppp)
}

check_lambda <- function (lambda, Nexpect) {
  N_will=sum(lambda$v*lambda$xstep*lambda$ystep)
  if(N_will<0.1){
    warning("the lambda is too small, it will take much time to generate a point pattern with at least one point")
  }
  if(N_will>(5*Nexpect)){
    maxlambda=5*Nexpect/(lambda$xstep*lambda$ystep)
    extri=which(lambda$v>maxlambda)
    lambda$v[extri]=maxlambda
    adn=Nexpect/sum(lambda$v*lambda$xstep*lambda$ystep)
    lambda$v=lambda$v*adn
  }
  return(lambda)
}