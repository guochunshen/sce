#' apply the same function on every member of a group
#' 
#' 
#' 
#' @param group factor, defination of unit of individuals used to fit model. e.g. species name, DBH class
#' 
#' 

applyGroup<-function(com,group,Fun,multicore=FALSE,mc.cores=1,...){
  
  if(is.null(group))
    group=as.factor(rep(1,com$N))
  grplevels=levels(group)
  
  for(i in 1:length(grplevels)){
    
    
  }
  
  
}