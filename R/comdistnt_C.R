#'
#'An C++ implementation of inter-community mean nearest taxon distance
#'
#'
#'


comdistnt_C<-function (comm, dis, abundance.weighted = FALSE, exclude.conspecifics = FALSE) 
{
  dat <- match.comm.dist(comm, dis)
  comm <- dat$comm
  dis <- dat$dist
  N <- dim(comm)[1]
  comm <- decostand(comm, method = "total", MARGIN = 1)
  if(abundance.weighted){
    x=t(apply(comm,1,function(x) x/sum(x)/2))
  }else{
    x=comm
    x[x!=0]=1
  }
  comdisnt=comdistntInner(N,dis,x,exclude.conspecifics)
  rownames(comdisnt) <- colnames(comdisnt) <- rownames(comm)
  return(as.dist(t(comdisnt)))
}