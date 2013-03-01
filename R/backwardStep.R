#' To automatically select a best fitted clustering model based on pvalue
#' 
#'@param fittedmodel a fm object represent a fitted clustering model for a point pattern
#'@param siglevel a marginal significance level used to decide which part of model is significance.
#'
#'@details
#' The \code{backwardStep} accept a \code{fm} object from \code{\link{fitCluster}}, then check whether
#' the trend of the model can be improved. Specifically, the improvement of trend is carried out by 
#' deleting insignificant environmental variables until all of then are significant or none of them left.
#' significant pvalues of each environmental variables were estimated by considering any possible significant clustering 
#' in the residual.
#' 
#' 
#' @seealso
#' \code{\link{fitCluster}}, \code{\link{updateCluster}}, \code{\link{varDecomp}}

backwardStep <- function(fittedmodel,siglevel=0.05){
  can_improve=TRUE
  while(can_improve){
    old_formula=as.character(attr(fittedmodel,"trend"))[2]
    fittedmodel=updateCluster(fittedmodel,siglevel=siglevel)
    new_formula=as.character(attr(fittedmodel,"trend"))[2]
    if(old_formula==new_formula)
      can_improve=FALSE
  }
  return(fittedmodel)
}