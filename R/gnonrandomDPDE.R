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
#'@examples
#'
#'r_samples=matrix(runif(100,0,4),nrow=50,ncol=2)
#'gnonrandomDPDE(r_samples,dtype="ptoe",k=1,area=100)
#'gnonrandomDPDE(r_samples,dtype="etoe",k=1,area=100)
#'


#TODO add the etoe estimator into the function
gnonrandomDPDE<-function(r_samples,dtype="ptoe",k,area,...){
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
    if(dtype=="ptoe"){
      ngpars=try(optim(c(n+1,2),fn=ngobf_ptoe,r=r,k=k,q=q,area=area,lower=c(n,1e-5),control=list(maxit=100000,fnscale=-1),...))
    }else{
      ngpars=try(optim(c(n+1,2),fn=ngobf_etoe,r=r,k=k,q=q,area=area,lower=c(n,1e-5),control=list(maxit=100000,fnscale=-1),...))
    }
    
    if(class(ngpars)=="try-error"){
      #browser()
      N=NA
    }else if(ngpars$convergence==0)
      N=ngpars$par[1]
    else
      N=NA
  }else{
    N=NA
  }
  return(N)
}

ngobf_ptoe=function(x,r,k,q,area){
  
  nind=x[1]
  lam=nind/area
  a=x[2]
  n=length(r)
  
  re=n*k*log(lam)+n*lgamma(k+a)+n*a*log(a)-n*lgamma(a)+
    sum(-(k+a)*log(pi*lam/q*r^2+a))
  return(re)
}


ngobf_etoe=function(x,r,k,q,area){
  nind=x[1]
  lam=nind/area
  a=x[2]
  n=length(r)
  
  re=n*k*log(lam)+n*lgamma(k+a+1)+n*(a+1)*log(a)-n*lgamma(a+1)+
    sum(-(k+a+1)*log(pi*lam/q*r^2+a))
  return(re)
}
