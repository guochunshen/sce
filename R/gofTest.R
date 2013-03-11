#'
#' a general goodness-of-fit (gof) test on the fitted cluster model
#'
#'@param a fittedmodel returned by \code{\link{fitCluster}} function
#'@param SFun a summary statistic of spaital point pattern
#'@param rRange the neighborhood distance range used to check performance of the model.
#'@param nsim number of simulated patterns.
#'@param r Numeric vector. The values of the argument r at which summary statistic should be evaluated.
#'@param ... other parameters passed to the summary statistic
#'
#'
#'@details
#'We could use a similar technique to do a goodness-of-fit test for the fitted log Gaussian Cox process. Consider e.g. the F-function. We can compute an estimate F_dat from the data and we can obtain estimated F-functions F_i, i=1,...,n from simulations of the fitted log Gaussian Cox process. We can then approximate the theoretical value of the F function by the mean F_mean of the simulated F_i, i=1,...,n. Our goodness of fit statistic then becomes
#'
#'D=integral from r_min to r_max of the squared distance (F_dat(r)-F_mean(r))^2.
#'
#'Similarly we obtain simulated values of D as
#'
#'D_i=integral from r_min to r_max of the squared distance (F_i(r)-F_mean(r))^2.
#'
#'Finally our goodness of fit p-value is  (m+1)/(n+1) where m is the number of D_i which are larger than D.
#'
#'The above procedure can be used with all the summary statistics. However the power of the goodness of fit test is weak with g since g was used for the model fitting (hence we should expect the theoretical g to be close to the g estimated from the data). Since K is just the integrated g, the goodness of fit test will also be weak with K. Nevertheless, we could carry out the goodness of fit test both with g, K, F, G, and J. 
#'
#'All the statistics are available in spatstat. Note however, that "inhomogeneous" versions of F, K, and G are not available so here we have to use the "stationary" versions.
#'
#'an insignificant pvalue means the fitted model has no significant different with the observed data.
#'
#'@return
#'a vector of pvalue, each pvalue corresponding to a edge correction method.
#'
#'
#'@examples
#'
#'#load one species data
#'data(testData)
#'sp1=subset(testData,testData$traits$species=="BEILPE")
#'
#'#fit a cluster model
#'fittedmodel=fitCluster(sp1,~elev+grad)
#'
#'#goodness-of-fit test
#'gofTest(fittedmodel,SFun=Fest,rRange=c(0,10),nsim=20,r=seq(0,20,1))
#'

gofTest<-function(fittedmodel,SFun,rRange=c(0,10),nsim=10,r=seq(0,20,1),...){
  #extract the real population data
  realdata=attr(fittedmodel,"data")
  
  #calculate the summary statistic of real data
  sm_real=SFun(realdata$com,r=r,...)
  sel=sm_real$r>=rRange[1] & sm_real$r<=rRange[2]
  sm_real=as.data.frame(sm_real)[sel,-c(1:2)]
  
  #calculate the summary statistic of simulated data
  sm_simu=list()
  for(i in 1:nsim){
    #generate a realization of a fitted model
    simudata=rCluster(fittedmodel,realdata$N,realdata$habitat)
    sm_simu[[i]]=as.data.frame(SFun(simudata,r=r,...))[sel,-c(1:2)]
    if(i==1)
      sm_mean=sm_simu[[i]]
    else
      sm_mean=sm_simu[[i]]+sm_mean
  }
  #calculate mean of the summary statistic
  sm_mean=sm_mean/nsim
  
  if(!is.null(dim(sm_real))){
    
    #calculate the D value of real data
    D_real=apply((sm_real-sm_mean)^2,2,sum2)
    
    #calculate the D value for each simulated data
    D_i=unlist(lapply(sm_simu,function(x){apply((x-sm_mean)^2,2,sum2)}))
    dim(D_i)=c(length(D_real),length(D_i)/length(D_real))
    
    #calculate pvalue
    m=apply(D_i>D_real,1,sum2)
  }else{
    #calculate the D value of real data
    D_real=sum2((sm_real-sm_mean)^2)
    
    #calculate the D value for each simulated data
    D_i=unlist(lapply(sm_simu,function(x){sum2((x-sm_mean)^2)}))
    
    #calculate pvalue
    m=sum2(D_i>D_real)
  }
  pvalues=(m+1)/(nsim+1)
  names(pvalues)=names(D_real)
  
  return(pvalues)
}

#a sum version of remove infinite value 
sum2<-function(x,na.rm=TRUE,inf.rm=TRUE){
  if(inf.rm){
    x=x[x!=Inf | x!=(-Inf)]
  }
  return(sum(x,na.rm=na.rm))
}