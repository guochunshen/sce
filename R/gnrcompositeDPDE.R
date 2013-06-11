#'
#' The general composite DPDE based on negative binomial distribution
#' 
#' 
#'@examples
#'ptoe_r=matrix(runif(100,1,4),nrow=50,ncol=2) 
#'etoe_r=matrix(runif(100,1,3),nrow=50,ncol=2)
#'
#'gnrcompositeDPDE(ptoe_r,etoe_r,k=2,area=100)
#'

gnrcompositeDPDE<-function(ptoe_r,etoe_r,k,area){
  N1=gnonrandomDPDE2(ptoe_r,dtype="ptoe",k,area)
  N2=gnonrandomDPDE2(etoe_r,dtype="etoe",k,area)
  N=sqrt(N1*N2)
  return(N)
}

gnrcompositeDPDE2<-function(ptoe_r,etoe_r,k,area){
  q=dim(ptoe_r)[2]
  
  #expected distance of event to event distance
  epr=mean(ptoe_r,na.rm=T)  
  eer=mean(etoe_r,na.rm=T)
  
  if(epr<eer){
    warning("invalid distance samples")
    return(NA)
  }
  
  #the aggregated paramter
  a=epr/(epr-eer)/2
  if(a<0.5){
    warning("invalid distance samples")
    return(NA)
  }
  
  #the lambda
  lam=q/pi*(sqrt(a)/epr * exp(lgamma(a-0.5)-lgamma(a)+lgamma(k+0.5)-lgamma(k)) )^2
  
  Nc=lam*area
  
  return(Nc)
}
