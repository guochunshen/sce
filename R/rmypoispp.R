#'
#' generate a poisson distribution point pattern, avoid 0 and very large number of point
#'
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