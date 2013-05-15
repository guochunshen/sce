#'
#'An C++ implementation of inter-community mean nearest taxon distance
#'
#'@param comm Community data matrix
#'@param interspecific distance matrix
#'@param abundance.weighted Should mean nearest taxon distances from each species to species in the other community be weighted by species abundance? (default = FALSE)
#'@param exclude.conspecifics Should conspecific taxa in different communities be exclude from MNTD calculations? (default = FALSE)
#'
#'@note it is an C++ implementation of the comdistnt function in picante package
#'
#'@examples
#'
#'library(testthat)
#'
#'re21=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = FALSE)
#'#test on the C++ version of the comdist function in the picante package
#'re22=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = FALSE)
#'#test the result
#'expect_equal(as.matrix(re21),as.matrix(re22))
#'re23=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = FALSE)
#'#test on the C++ version of the comdist function in the picante package
#'re24=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = FALSE)
#'expect_equal(as.matrix(re13),as.matrix(re14))
#'
#'re21=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = TRUE)
#'#test on the C++ version of the comdist function in the picante package
#'re22=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = TRUE)
#'#test the result
#'expect_equal(as.matrix(re21),as.matrix(re22))
#'re23=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = TRUE)
#'#test on the C++ version of the comdist function in the picante package
#'re24=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = TRUE)
#'expect_equal(as.matrix(re13),as.matrix(re14))
#'
#'
#'


comdistnt_C<-function (comm, dis, abundance.weighted = FALSE, exclude_conspecifics = FALSE, cal_pairs=NULL) 
{
  dat <- match.comm.dist(comm, dis)
  comm <- dat$comm
  dis <- dat$dist
  N <- dim(comm)[1]
  S = dim(comm)[2]
  comm <- decostand(comm, method = "total", MARGIN = 1)
  if(abundance.weighted){
    x=t(apply(comm,1,function(x) x/sum(x)/2))
  }else{
    x=comm
    x[x!=0]=1
    x=x/apply(x,1,sum)/2
  }
  if(is.null(cal_pairs)){
    cal_pairs=matrix(TRUE,nrow=N,ncol=N)
  }
  if(any(is.nan(x))){
    nax=apply(x,1,function(y) any(is.nan(y)))
    cal_pairs[(rep(nax,each=N) | rep(nax,times=N))]=FALSE
    x[is.nan(x)]=0
  }
  #comdisnt=comdistntInner(N,dis,x,exclude_conspecifics,cal_pairs)
  comdisnt=rep(0,N^2)
  comdisnt=.C("comdistnt",as.integer(N),as.integer(S),as.double(dis),as.double(x),as.integer(cal_pairs),
              as.integer(exclude_conspecifics),as.double(comdisnt))[[7]]
  dim(comdisnt)=c(N,N)
  #set the uninteresting community pair equals to NA
  comdisnt[!cal_pairs]=NA
  
  rownames(comdisnt) <- colnames(comdisnt) <- rownames(comm)
  re=as.dist(comdisnt)
  re[is.na(re)]=NA
  return(re)
}