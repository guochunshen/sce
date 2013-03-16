#' To test whether there is any aggregative residual existed in the unexplianed pattern by habitat
#' 
#' @param fittedmodel a \link{fm} object represented a fitted cluster point process model for a point pattern.
#' @param nsim number of relazations used to simulate point pattern accroding to the model.
#' @param r distance range used in goodness of fit test.
#' @param edcor edge correction method. see details in \link{pcf}
#'
#'@details
#' This is a goodness-of-fit test based on the heterogeneous pari correlation function.
#'
#'@examples
#'data(testData)
#'
#' sp1=subset(testData,testData$traits$species=="ACALDI")
#' 
#' fm=fitCluster(sp1,~elev+grad)
#' 
#' sigAggreResidualTest(fm,nsim=10,r=0:30)
#'

sigAggreResidualTest=function(fittedmodel,nsim=10,r=seq(0,60,2),edcor="translate"){
  data=attr(fittedmodel,"data")
  xy=data.frame(x=data$com$x,y=data$com$y)
  data.ppm=attr(fittedmodel,"fittedmodel")
  bw=attr(fittedmodel,"ctlpars")$bw
  npp=nsim+1
  
  #since the translate edge correction need very large memory in its calculation, thus using border edge correction
  #instand is the number of points large than 3000
#   if(data$N<3000){
#     e=envelope(data.ppm,fun=pcf,nsim=nsim,verbose=FALSE,savefuns=TRUE,correction=edcor,r=r)
#   }else{
#     e=envelope(data.ppm,fun=pcf,nsim=nsim,verbose=FALSE,savefuns=TRUE,correction="border",r=r)
#   }
#   Kfuns=attr(e,"simfuns")
#   n=length(e$obs)
#   Kfuns[,1]=e$obs
  
  lambda_loc=predict.ppm(data.ppm,locations=xy)
  pcf_obs=pcf_adaptive(data$com,lambda=lambda_loc,maxt=max(r),bw=bw,adaptive =FALSE)
  
  pcf_r=pcf_obs$r
  pcf_obs=pcf_obs$trans
  #plot(x=pcf_r,y=pcf_obs)
  
  Kfuns=matrix(nrow=length(pcf_obs),ncol=npp)
  lambda=predict.ppm(data.ppm)
  for(i in 1:nsim){
    data_simu=rpoispp(lambda)
    Kfuns[,i+1]=pcf_adaptive(data_simu,lambda=lambda_loc,maxt=max(r),bw=bw,adaptive =FALSE)$trans
    while(any(is.infinite(Kfuns[,i+1]))){
      data_simu=rpoispp(lambda)
      Kfuns[,i+1]=pcf_adaptive(data_simu,lambda=lambda_loc,maxt=max(r),bw=bw,adaptive =FALSE)$trans
    }
     
    #lines(x=pcf_r,y=Kfuns[,i+1],col=2)
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
