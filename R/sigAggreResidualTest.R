#' To test whether there is any aggregative residual existed in the unexplianed pattern by habitat
#' 
#' @param fittedmodel a \link{fm} object represented a fitted cluster point process model for a point pattern.
#' @param nsim number of relazations used to simulate point pattern accroding to the model.
#' @param r distance range used in goodness of fit test.
#' @param edcor edge correction method. see details in \link{pcf}
#'
#'@details
#' This is a goodness-of-fit test based on the heterogeneous pari correlation function. 
#' since \code{\link{pcf_adaptive}} is good at lower memory and worse on scater data point pattern 
#' (e.g. number of individuals below than 50). Thus we use \code{\link{pcf_adaptive}} on abundant species,
#' and use \code{\link{pcf}} on rare species.
#'
#'@examples
#'data(testData)
#'
#' sp1=subset(testData,testData$traits$species=="ACALDI")
#' 
#' fm=fitCluster(sp1,~elev+grad)
#' 
#' #significant residual
#' sigAggreResidualTest(fm,nsim=10,r=0:30)
#' 
#' sp2=subset(testData,testData$traits$species=="HIRTTR")
#'
#' fm2=fitCluster(sp2,~elev+grad)
#' 
#' sigAggreResidualTest(fm2,nsim=10,r=0:30)
#' 

sigAggreResidualTest=function(fittedmodel,nsim=10,r=seq(0,60,2),edcor="translate",ntry=10){
  data=attr(fittedmodel,"data")
  xy=data.frame(x=data$com$x,y=data$com$y)
  data.ppm=attr(fittedmodel,"fittedmodel")
  bw=attr(fittedmodel,"ctlpars")$bw
  npp=nsim+1
  lambda=predict.ppm(data.ppm)
  #if we get extrem value in the lambda, it means the ppm model is not well fitted
  #thus we can adjust the extrem values in the lambda to avoid memory problem
  lambda=check_lambda(lambda,data$N)
  
  pcf_obs=get_pcf(data$com,lambda,xy,r,edcor,bw,data$N) 
  
  pcf_r=pcf_obs$r
  pcf_obs=pcf_obs$trans
  #plot(x=pcf_r,y=pcf_obs)
  
  Kfuns=matrix(nrow=length(pcf_obs),ncol=npp)
  
  
  for(i in 1:nsim){
    data_simu=rmypoispp(lambda,data$N)
    Kfuns[,i+1]=get_pcf(data_simu,lambda,xy,r,edcor,bw,data$N)$trans
    localn=1
    while(any(is.infinite(Kfuns[,i+1])) & localn<ntry){
      data_simu=rmypoispp(lambda,data$N)
      Kfuns[,i+1]=get_pcf(data_simu,lambda,xy,r,edcor,bw,data$N)$trans
      localn=localn+1
    }
    
    #lines(x=pcf_r,y=Kfuns[,i+1],col=2)
  }
  #if there is always infinite value contained in the pcf, just return a pvalue equals to 1
  if(localn>=ntry){
    pvalue=1
    names(pvalue)="aggreRes"
    return(pvalue)
  }
  Kfuns[,1]=pcf_obs
  n=length(pcf_obs)
  
  
  hsum=apply(Kfuns,1,sum,na.rm=TRUE)
  hmean=(hsum-Kfuns)/nsim
  
  ui=apply((Kfuns-hmean)^2,2,sum,na.rm=TRUE)
  pvalue=1-sum(ui[1]>ui[-1])/nsim
  
  
  #   
  #   
  #   ui=vector("numeric",npp)
  #   for (k in 2:(n-1)){
  #     Kit=Kfuns[k,]
  #     t.rslt=step.ui(Kit,2,npp)
  #     ui = ui+t.rslt
  #   }
  #   pvalue=calc.pval(ui)
  names(pvalue)="aggreRes"
  return(pvalue)
}


get_pcf <- function (com,lambda,xy,r,edcor,bw,N) {
  lambda_loc=interp.im(lambda,x=com$x,y=com$y)
  if(N<500){
    pcf_obs=try(pcf(com,lambda=lambda_loc,r=seq(0,max(r),length.out=101),correction=edcor))
    pcf_obs=pcf_obs[-1,]
  }else{
    pcf_obs=pcf_adaptive(com,lambda=lambda_loc,maxt=max(r),bw=bw,adaptive =FALSE)
  }
  return(pcf_obs)
}


step.ui <- function(Kit, delta.t, npp, c=0.25)
{
  Ki.bar <- (sum(Kit) - Kit)/(npp-1)
  
  this.ui <- (Kit^c-Ki.bar^c) *
    (Kit^c-Ki.bar^c) * delta.t
  
  
  ui.step <- vector("numeric", npp)
  ui.step <- ui.step + this.ui
  
  return(ui.step)
}

calc.pval <- function(ui)
{
  # fetch the ui for the observed pattern
  my.ui <- ui[1]
  # calculate the number of values greater than the observed
  my.pval <- (length(ui[ui>=my.ui]))/(length(ui))
  
  return(my.pval)
}
