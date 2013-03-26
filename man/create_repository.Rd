\name{create_repository}
\alias{create_repository}
\title{a result repository interface, it designed for long time calculation and many replications.}
\usage{
  create_repository(identity,schedule,repository_path=".")
    read_repository(identity, repository_path=".")
    getUnfinishedJobIds(rep) getUnfinishedJobNames(rep)
    writeResultById(rep,jobid,result)
    remove_repository(rep,repository_path=".")
}
\arguments{
  \item{identity}{a character represent the uniqueness of
  the repository object (the file name of an RData file).}

  \item{schedule}{a vector represent all steps of the
  calculation}

  \item{repository_path}{the location of the repository
  file}

  \item{rep}{an repository object}

  \item{jobid}{a numeric id represent the location of the
  current progress}
}
\description{
  a result repository interface, it designed for long time
  calculation and many replications.
}
\details{
  every job was stored seperately, and managed by
  repository object
}
\examples{
#create a reposite at current working director
rp1=create_repository("test",paste("job",1:5),".")

#get unfinished job ids
getUnfinishedJobIds(rp1)


#write one result into the repository
writeResultById(rp1,2,data.frame(1))


#get unfinished job names
getUnfinishedJobNames(rp1)

#read repository into other object
rp2=read_repository("test")

#you can get the same unfinished job names
getUnfinishedJobNames(rp2)

#get all results in the repository
getAllResults(rp1)

#finally, you can delete the repository
remove_repository(rp2)
}

