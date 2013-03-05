#' To update a model manually or automatically
#'
#' @param fittedmodel a fm object represent a fitted clustering model for a point pattern
#' @param trend a \code{\link{formula}} object represent the structure of the intensity trend.
#' if it is null, the most insignificant habitat will be removed from the formula.
#' @param siglevel a marginal significance level used to decide which part of model is significance.
#' 
#' @seealso \code{\link{fitCluster}}, \code{\link{sigHabitatTest}}, \code{\link{sigAggreResidualTest}}
#' 
#'

updateCluster <- function(fittedmodel,trend=NULL,siglevel=0.05){
  
  pvalues=attr(fittedmodel,"pvalues")
  if(is.null(pvalues)){
    pvalues=sigTestofCluster(fittedmodel)
    attr(fittedmodel,"pvalues")=pvalues
  }
  
  if(class(pvalues)=="try-error"){
    warning("can't evaluate the model caused by error in calculation of pvalues")
    return(fittedmodel)
  }
  
  modeltype=attr(fittedmodel,"modeltype")
  data=attr(fittedmodel,"data")
  ctlpars=attr(fittedmodel,"ctlpars")
  cluster=attr(fittedmodel,"modeltype")
  
  #1. whether the model can be improve?
  #1.1 delete any insignificance habitat first
  hasinsig=length(pvalues)>1 & any(pvalues[-1]>siglevel)
  
  #2. if yes, delete some part of the model and given a new formular
  if(hasinsig & is.null(trend)){
    del=which(pvalues[-1]==max(pvalues[-1],na.rm=T))[1]
    allterms=names(pvalues[-1])[-del]
    if(length(allterms)>0)
    trend=as.formula(paste("~ ",paste(allterms,collapse="+")))
    else
      trend=as.formula("~1")
    #3. update the fitted model based the new formular
    fittedmodel=fitCluster(data,trend,cluster,sigTest=TRUE,ctlpars=ctlpars)
  }
  
  return(fittedmodel)
}