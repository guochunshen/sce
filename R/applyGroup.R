#' apply the same function on every member of a group
#' 
#' 
#' @param com a \code{\link{scp}} object.
#' @param group factor, defination of unit of individuals used to fit model. e.g. species name, DBH class
#' @param Fun a function name that will apply on every element of a group
#' @param verbpro logical flage to show the verbose of process of calculation
#' @param ... other parameters passed to the \code{Fun} function
#' 
#' 
#' @return a list of results that returned by the function
#' 


applyGroup<-function(com,group,Fun,verbpro=TRUE,...){
  
  if(is.null(group))
    stop("please specify a valid group index")
  if(length(group)!=com$N)
    stop("Length of group and total community size are not equal")
  
  grplevels=unique(group)
  re=list()
  for(i in 1:length(grplevels)){
    if(verbpro){
      print(paste("starting the", i, "th element in the group"))
    }
    sel=grplevels[i]==group
    subdata=subset(com,sel)
    
    re[[i]]=Fun(subdata, ...)
  }
  return(re)
}
