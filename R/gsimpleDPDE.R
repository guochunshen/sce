#'
#' The general simple distance based population density estimator (DPDE)
#'
#'@param r_samples a matrix of distance, distances in each column comes from one section of q equal-angle sectors.
#'@param k the order of nearest neighbor
#'@param area the area of region that want to estimate number of individuals
#'
#'
#'@examples
#'
#'r_samples=matrix(runif(100,0,4),nrow=50,ncol=2)
#'gsimpleDPDE(r_samples,k=1,area=100)
#'
#'


gsimpleDPDE<-function(r_samples,k,area){
  q=dim(r_samples)[2]
  del=apply(r_samples,1,function(x) any(is.na(x)))
  r_samples=r_samples[!del,]
  n=dim(r_samples)[1]
  if(is.null(n))
    n=length(r_samples)
  N1=area/sum(r_samples^2)*q*(q*n*k-1)/pi
  return(N1)
}