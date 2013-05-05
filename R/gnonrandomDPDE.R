#'
#'The general nonrandom distance based population density estimator (DPDE)
#'
#'@param r_samples neighborhood distance matrix. The number of column represents the number of equal angle sectors used in distance sampling
#'@param dtype type of neighborhood distance in \code{r_sample}. it could be "ptoe" (point to event distance) or "etoe" (event to event distance).
#'@param k the order of neighborhood distance
#'@param area area the area in which population size want to be estimated
#'
#'@details
#'This population density esimator is based on a assumption of population density distribution across quadrat following a negative
#'binomial distribution. In this case, the forms of estimator based on point to event distances and event to event distances are
#'different.
#'
#'


#TODO add the etoe estimator into the function
gnonrandomDPDE<-function(r_samples,dtype="ptoe",k,area){
  q=dim(r_samples)[2]
  del=apply(r_samples,1,function(x) any(is.na(x)))
  r_samples=r_samples[!del,]
  n=dim(r_samples)[1]
  if(is.null(n))
    n=length(r_samples)
  r=as.vector(r_samples)
  #remove 0 neighborhood distance
  r=r[r!=0]
  if(length(r)>0){
    ngpars=try(optim(c(0.03,2),fn=ngobf,r=r,k=k,q=q,control=list(maxit=100000,fnscale=-1)))
    if(class(ngpars)=="try-error"){
      #browser()
      N=NA
    }else if(ngpars$convergence==0)
      N=ngpars$par[1]*area
    else
      N=NA
  }else{
    N=NA
  }
  return(N)
}

ngobf=function(x,r,k,q){
  #x=c("lambda","a")
  if(any(x<0))
    return(-Inf)
  lam=x[1]
  a=x[2]
  n=length(r)
  #re=n*log(lam)+sum(log(r))-(k+1)*sum(log(1+(pi*lam*r^2)/k))
  re=n*log(lam)+n*lgamma(k+a)+n*a*log(a)-n*lgamma(a)+
    sum(log(r)+(k-1)*log(pi*lam/q*r^2)-(k+a)*log(pi*lam/q*r^2+a))
  return(re)
}
