#' To automatically select a best fitted clustering model based on pvalue
#' 
#'@param fittedmodel a fm object represent a fitted clustering model for a point pattern
#'@param siglevel a marginal significance level used to decide which part of model is significance.
#'
#'

backwardStep <- function(fittedmodel,siglevel=0.05){
  can_improve=TRUE
  while(can_improve){
    old_formula=attr(fittedmodel,"trend")
    fittedmodel=updateCluster(fittedmodel,siglevel)
    new_formula=attr(fittedmodel,"trend")
    if(old_formula==new_formula)
      can_improve=FALSE
  }
  return(fittedmodel)
}