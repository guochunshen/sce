#'To estimate the proportion of variance explained by habitat, internal cluster and random process from a given fitted cluster model
#'
#'@param fittedmodel a \code{\link{fm}} object representing a fitted clustering model
#'@param r spatial range used in pair correlation function to caluclate the proportaion of variance
#'@param R the quadrat scale used to evaluate the variance of poisson noise
#'@param delta use delta x delta grid for numerical integration.
#'@param simple if it is true, only PVH and variance of habitat and spatial autocorrelation were calculated
#'
#'@return a vector of length 6 contains \code{"PVH","PVR","varHabitat","varInterCluster","varPoi","varNonrand"}. Specifically,
#'  "PVH": proportion of variance explained by Habitat, thus proportion of variance explained
#'  by internal clustering processes is \bold{PVI}=1-PVH
#'  "PVR": proportion of variance explained by poisson Random process, 
#'  thus proportion of variance explained by the nonrandom processes (\bold{PVNR}) is 1-PVR
#'  "varHabitat" and "varInterCluster" are the variance explained by habtiat and internal cluster
#'  process in the log intensity function of the point pattern; "varPoi" and "varNonrand" are the variances
#'  explained by poisson random and nonrandom processes in the intensity function of the point pattern.
#'


varDecomp<-function(fittedmodel,r=c(0:80),R=10,delta=1,simple=FALSE){
  data.ppm=attr(fittedmodel,"fittedmodel")
  trendmap=predict(data.ppm, type="trend")
  varHabitat=var(log(as.numeric(trendmap$v)))
  varInterCluster=fittedmodel[2]
  PVH=varHabitat/(varHabitat+varInterCluster)
  
  if(simple){
    re=c(PVH,NA,varHabitat,varInterCluster,NA,NA)
    names(re)=c("PVH","PVR","varHabitat","varInterCluster","varPoi","varNonrand")
    return(re)
  }
  
  data.ppp=attr(fittedmodel,"data")$com
  pcfvalues=pcf(data.ppp,r=r)
  #pcf may return infinite value at zero. We replace Inf by the maximal value less than Inf.
  pcfvalues$trans[pcfvalues$trans==Inf]=max(pcfvalues$trans[pcfvalues$trans!=Inf])
  #pcf is close to one at maximum distance meters. We want to not overestimate the
  #proportion of variance due to habitat or seed dispersal so I insert the value of one at the maximal value of r.
  pcfvalues$trans[length(pcfvalues$r)]=1
  
  rho=data.ppp$n/area.owin(data.ppp$win)
  
  #R=10, consider 10 by 10 m square
  #delta=1, we use 1 x 1 grid for numerical integration
  
  #consider variance of total number of trees in a r x r square.
  
  #Poisson variance:
  VarPoi=rho*R^2
  
  #Variance due to seed dispersal and environment is computed using numerical quadrature.
  X=ppp(x=c(),y=c(),window=owin(xrange=c(0,R),yrange=c(0,R)))#this is just an empty dummy point pattern
  Q=quadscheme(X,nd=floor(R/delta))#this function creates the numerical integration quadrature points and weights
  
  rr <- pairdist(union.quad(Q), squared = FALSE)#computes all pairwise distances between quadrature points
  integrand=pcf.lookup(rr,pcfvalues,length(pcfvalues$r),1)-1#evaluate integrand pcf()-1 for all pairs of quadrature points
  integrand=matrix(integrand,ncol=dim(rr)[1])#convert into matrix
  doubleintegral=Q$w%*%integrand%*%Q$w#compute double integral as matrix product.
  #variance due to habitat and seed dispersal (spatially structured variation)
  VarNonrand=rho^2*doubleintegral
  PVR=VarPoi/(VarPoi+VarNonrand)
  
  re=c(PVH,PVR,varHabitat,varInterCluster,VarPoi,VarNonrand)
  names(re)=c("PVH","PVR","varHabitat","varInterCluster","varPoi","varNonrand")
  return(re)
}


pcf.lookup=function(r,fittedpcf,maxindex,deltar){
  #maxindex is the maximal index for fittedpcf$r. deltar is the distance between  #values in fittedpcf$r
  index=pmin(ceiling(r/deltar),maxindex-1)+1
  
  return(fittedpcf$trans[index])
}
