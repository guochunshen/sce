#'
#' The general composite distance based population density estimator (DPDE)
#'
#'@param ptoe_r neighborhood distances matrix from the focal point/location to its nearest neighbor, the number of column represent the number of equal angle sector used in sampling distances
#'@param etoe_r neighborhood distances matrix from the focal event to its nearest neighbor
#'@param k the order of neighbor
#'@param area the area in which population size want to be estimated
#'
#'@examples
#'ptoe_r=matrix(runif(100,0,4),nrow=50,ncol=2)
#'#ptoe_r and etoe_r can have different sample sizes
#'etoe_r=matrix(runif(120,0,4),nrow=60,ncol=2)
#'
#'gcompositeDPDE(ptoe_r,etoe_r,k=2,area=100)
#'

gcompositeDPDE<-function(ptoe_r,etoe_r,k,area){
  q=dim(ptoe_r)[2]
  del=apply(ptoe_r,1,function(x) any(is.na(x)))
  ptoe_r=ptoe_r[!del,]
  del=apply(etoe_r,1,function(x) any(is.na(x)))
  etoe_r=etoe_r[!del,]
  npe=dim(ptoe_r)[1]
  nee=dim(etoe_r)[1]
  if(is.null(npe)){
    npe=length(ptoe_r)
    nee=length(etoe_r)
  }
  Nc=q*exp(lgamma(npe*q*k)+lgamma(nee*q*k)-lgamma(npe*q*k-1/2)-lgamma(nee*q*k-1/2))/pi*area/sqrt(sum(ptoe_r^2))/sqrt(sum(etoe_r^2))
  return(Nc)
}
