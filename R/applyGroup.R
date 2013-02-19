#' apply the same function on every member of a group
#' 
#' 
#' @param com a \code{\link{scp}} object.
#' @param group factor, defination of unit of individuals used to fit model. e.g. species name, DBH class
#' 
#' 

applyGroup<-function(com,group,Fun,multicore=FALSE,mc.cores=1,...){
  
  if(is.null(group))
    stop("please specify a valid group index")
  
  grplevels=unique(groups)
  for(i in 1:length(grplevels)){
    
  }
  
  
}