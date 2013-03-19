#' a result repository interface, it designed for long time calculation and many replications. 
#'
#' @usage
#' create_repository(identity,schedule,repository_path=".")
#' read_repository(identity, repository_path=".")
#' getUnfinishedJobIds(rep)
#' getUnfinishedJobNames(rep)
#' writeResultById(rep,jobid,result)
#' remove_repository(rep,repository_path=".")
#'
#'@param identity a character represent the uniqueness of the repository object (the file name of an RData file). 
#'@param schedule a vector represent all steps of the calculation
#'@param repository_path the location of the repository file
#'@param rep an repository object
#'@param jobid a numeric id represent the location of the current progress
#'
#'@details
#' every job was stored seperately, and managed by repository object
#'
#'@examples
#'
#'#create a reposite at current working director
#'rp1=create_repository("test",paste("job",1:5),".")
#'
#'#get unfinished job ids
#'getUnfinishedJobIds(rp1)
#'
#'
#'#write one result into the repository
#'writeResultById(rp1,2,data.frame(1))
#'
#'
#'#get unfinished job names
#'getUnfinishedJobNames(rp1)
#'
#'#read repository into other object
#'rp2=read_repository("test") 
#'
#'#you can get the same unfinished job names
#'getUnfinishedJobNames(rp2)
#'
#'#finally, you can delete the repository 
#'remove_repository(rp2)
#'

create_repository<-function(identity,schedule,repository_path="."){
  #create the object
  n=length(schedule)
  rep=data.frame(jobnames=schedule,status=rep(FALSE,n))
  attr(rep,"identity")=identity
  #full path location of the file and the name of the file
  attr(rep,"filepath")=file.path(repository_path,paste(identity,".RData",sep=""))
  attr(rep,"njobs")=n
  class(rep)=c("repository","data.frame")
  
  #save the object into a seperate file
  save(rep,file=attr(rep,"filepath"))
  return(rep)
}

read_repository<-function(identity, repository_path=".",recheck=TRUE){
  file=file.path(repository_path,paste(identity,".RData",sep=""))
  if(file.exists(file)){
    load(file)
    if(recheck){
      n=attr(rep,"njobs")
      basepath=strsplit(attr(rep,"filepath"),".RData")[[1]]
      #check individual job files in the repository
      for(i in 1:n){
       jobfile=paste(basepath,"_",rep$jobnames[i],".RData")
       rep$status[i]=file.exists(jobfile)
      }
    }
    
    return(rep)
  }else{
    return(FALSE)
  }
}
# 
# update_repository<-function(rep){
#   if(!inherits(rep,"repository"))
#     stop("rep should be a repository object")
#   save(rep,file=attr(rep,"filepath"))
# }


#return the schedule numeric ids of unfinished jobs
getUnfinishedJobIds<-function(rep){
  ids=which(!rep$status)
  return(ids)
}

#return the schedule character name of unfinished jobs
getUnfinishedJobNames<-function(rep){
  ids=getUnfinishedJobIds(rep)
  jobnames=rep$jobnames[ids]
  return(jobnames)
}


writeResultById<-function(rep,jobid,result){
  if(!inherits(rep,"repository"))
    stop("rep should be a repository object")
  save(result,file=getJobfilename(rep,jobid))
  return(NULL)
}

readAllResults<-function(rep){
  n=attr(rep,"njobs")
  re=list()
  for(i in 1:n){
    jobfile=getJobfilename(rep,i)
    if(file.exists(jobfile)){
      load(jobfile)
      re[[i]]=result
    }else{
      warning(paste("the",i," th job has not been done"))
    }
  }
  return(re)
}

removeEmptyResult<-function(rep){
  n=attr(rep,"njobs")
  re=list()
  for(i in 1:n){
    jobfile=getJobfilename(rep,i)
    if(file.exists(jobfile)){
      load(jobfile)
      if(is.null(result)){
        #file.remove(jobfile)
        print(paste("the result of",rep$jobnames[i],"job is empty, removed!"))
      }
    }
  }
}

# 
# pushResultByName<-function(rep,progressName,result){
#   if(!inherits(rep,"repository"))
#     stop("rep should be a repository object")
#   rep[progressName]=result
#   return(rep)
# }

remove_repository<-function(rep,repository_path="."){
    file.remove(attr(rep,"filepath"))
    n=attr(rep,"njobs")
    for(i in 1:n){
      jobfile=getJobfilename(rep,i)
      if(file.exists(jobfile))
        file.remove(jobfile)
    }
}

print.repository<-function(x){
  cat("the",attr(x,"identity"),"repository.\n")
  cat(dim(x)[1],"jobs in total. \n")
  cat(sum(!x$status),"jobs are unfinished.")
}

getJobfilename<-function(rep,jobid){
  basepath=strsplit(attr(rep,"filepath"),".RData")[[1]]
  jobfile=paste(basepath,"_",rep$jobnames[jobid],".RData")
  return(jobfile)
}
