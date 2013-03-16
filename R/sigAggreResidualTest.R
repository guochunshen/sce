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
  data.ppm=attr(fittedmodel,"fittedmodel")
  #since the translate edge correction need very large memory in its calculation, thus using border edge correction
  #instand is the number of points large than 3000
  if(data$N<3000){
    e=envelope(data.ppm,nsim=nsim,verbose=FALSE,savefuns=TRUE,correction=edcor,r=r)
  }else{
    e=envelope(data.ppm,nsim=nsim,verbose=FALSE,savefuns=TRUE,correction="border",r=r)
  }
  
  npp=nsim+1
  Kfuns=attr(e,"simfuns")
  n=length(e$obs)
  Kfuns[,1]=e$obs
  ui=vector("numeric",npp)
  for (k in 2:(n-1)){
    Kit=Kfuns[k,]
    t.rslt=step.ui(Kit,2,npp)
    ui = ui+t.rslt
  }
  pvalue=calc.pval(ui)
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
