#' apply the same function on every member of a group, multicore is supported on Linux with multicores CPU
#' 
#' 
#' @param com a \code{\link{scp}} object.
#' @param group factor, defination of unit of individuals used to fit model. e.g. species name, DBH class
#' @param Fun a function name that will apply on every element of a group
#' @param verbpro logical flage to show the verbose of process of calculation
#' @param mc.cores number of cores used in parallel calculation on Linux
#' @param ... other parameters passed to the \code{Fun} function
#' @param repository an repository object that can save middle result into the repository file. 
#'        see \code{\link{create_repository}} for detail. it only usefull to run long term calculation 
#'        if the calculation environmental is not stable and will crash sometime.
#' 
#' @return a list of results that returned by the function
#' 


applyGroup<-function(com,group,Fun,verbpro=FALSE,mc.cores=5,multicore=FALSE,repository=NULL,...){
  
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
  if(is.null(repository)){
    joblist=1:nls
  }else{
    print("continuing doing unfinished jobs in the given repository")
    joblist=getUnfinishedJobIds(repository)
  }
  
  if(length(joblist)>0){
    if(multicore){
      re=mclapply(joblist,applyOneGroup,grplevels=grplevels,verbpro=verbpro,group=group,Fun=Fun,com=com,
                  repository=repository,mc.cores=mc.cores,mc.preschedule = FALSE,...)  
    }
    else
      re=lapply(joblist,applyOneGroup,grplevels=grplevels,verbpro=verbpro,group=group,Fun=Fun,com=com,repository=repository,...)
  }else if(!is.null(repository)){
    re=readAllResults(repository)
  }
  
  names(re)=as.character(grplevels)
  return(re)
}

applyOneGroup <- function ( i, grplevels, verbpro, group, Fun, com, repository, ...) {
    if(verbpro){
      print(paste("starting the", i, "th element in the group"))
    }
    sel=grplevels[i]==group
    subdata=subset(com,sel)
    
    re=Fun(subdata, ...)
    if(verbpro){
      print(paste("the", i, "th element finished"))
    }
    if(!is.null(repository)){
      repository=writeResultById(repository,i,re)
    }
    return(re)
}



