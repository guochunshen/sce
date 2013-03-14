#' apply the same function on every member of a group, multicore is supported on Linux with multicores CPU
#' 
#' 
#' @param com a \code{\link{scp}} object.
#' @param group factor, defination of unit of individuals used to fit model. e.g. species name, DBH class
#' @param Fun a function name that will apply on every element of a group
#' @param verbpro logical flage to show the verbose of process of calculation
#' @param mc.cores number of cores used in parallel calculation on Linux
#' @param ... other parameters passed to the \code{Fun} function
#' 
#' 
#' @return a list of results that returned by the function
#' 


applyGroup<-function(com,group,Fun,verbpro=FALSE,mc.cores=5,multicore=FALSE,...){
  
  if(is.null(group))
    stop("please specify a valid group index")
  if(length(group)!=com$N)
    stop("Length of group and total community size are not equal")
  
  grplevels=unique(group)
  nls=nlevels(group)
  if(multicore & Sys.info()["sysname"]=="Windows" ){
   warning("Current version of R only support parallel computation under linux, thus multicore is disabled")
   multicore=FALSE
  }
  
  if(multicore){
    re=mclapply(1:nls,applyOneGroup,grplevels=grplevels,verbpro=verbpro,group=group,Fun=Fun,com=com,
                mc.cores=mc.cores,...)  
  }
  else
    re=lapply(1:nls,applyOneGroup,grplevels=grplevels,verbpro=verbpro,group=group,Fun=Fun,com=com,...)
  
  names(re)=as.character(grplevels)
  return(re)
}

applyOneGroup <- function ( i, grplevels, verbpro, group, Fun, com, ...) {
    if(verbpro){
      print(paste("starting the", i, "th element in the group"))
    }
    sel=grplevels[i]==group
    subdata=subset(com,sel)
    
    re=Fun(subdata, ...)
    if(verbpro){
      print(paste("the", i, "th element finished"))
    }
    return(re)
}



