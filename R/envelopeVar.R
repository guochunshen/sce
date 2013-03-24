#'
#'calculate envelope of variances explained different fractions of the model
#'
#'@param fittedmodel a fm object returned by the \code{\link{fittedmodel}} function
#'@param nsim number of simulations used to generate the envelope
#'@param conf_level significance level of the confidence interval
#'@param r,R,delta parameters used in the \code{\link{varDecomp}} function
#'@param simple a logical flag in \code{\link{varDecomp}} to calculate the PVH only if it is true
#'

envelopeVar<-function(fittedmodel,nsim=9,conf_level=0.95,r=c(0:80),R=10,delta=1,simple=FALSE){
  #extract the real population data
  realdata=attr(fittedmodel,"data")
  
  #extract parameters
  ctlpars=attr(fittedmodel,"ctlpars")
  trend=attr(fittedmodel,"trend")
  #calculate the summary statistic of simulated data
  sm_simu=list()
  for(i in 1:nsim){
    #generate a realization (an scp object) of a fitted model
    simudata=rCluster(fittedmodel,realdata$N,realdata$habitat,ntry=3)
    #if there is some error appeared in the point pattern simulation, just return a NULL pvalues
    if(inherits(simudata,"try-error")){
        return(NULL)
    }
    sm_simu[[i]]=varDecomp(fitCluster(simudata,trend,sigTest=FALSE,ctlpars=ctlpars),r=r,R=R,delta=delta,simple=simple)
  }
  varnames=names(sm_simu[[1]])
  
  sm_simu=unlist(sm_simu)
  dim(sm_simu)=c(length(varnames),nsim)
  deli=apply(sm_simu,1,function(x) all(is.na(x)))
  sm_simu=sm_simu[!deli,]
  varnames=varnames[!deli]
  #calculate confidence of the variance
  re=apply(sm_simu,1,confij,conf_level=conf_level)
  #calculate the stand divation
  re_sd=apply(sm_simu,1,sd,na.rm=T)
  re=rbind(re,re_sd)
  colnames(re)=varnames
  return(re)
}

confij=function(x, conf_level=0.95){
  x=sort(x)
  n=length(x)
  bottomi=floor(n*(1-conf_level))
  if(bottomi==0)
    bottomi=1
  
  upi=ceiling(n*conf_level)
  if(upi>n)
    upi=n

  return(c(x[bottomi],x[upi]))
}