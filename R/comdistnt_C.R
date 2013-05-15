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
#'data(phylocom)
#'re11=comdist(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
#'#test on the C++ version of the comdist function in the picante package
#'re12=comdist_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
#'#test the result
#'expect_equal(re11,re12)
#'re13=comdist(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE)
#'#test on the C++ version of the comdist function in the picante package
#'re14=comdist_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE)
#'expect_equal(re13,re14)
#'
#'re21=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = FALSE)
#'#test on the C++ version of the comdist function in the picante package
#'re22=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude_conspecifics = FALSE)
#'#test the result
#'expect_equal(as.matrix(re21),as.matrix(re22))
#'re23=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = FALSE)
#'#test on the C++ version of the comdist function in the picante package
#'re24=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude_conspecifics = FALSE)
#'expect_equal(as.matrix(re13),as.matrix(re14))
#'
#'re21=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = TRUE)
#'#test on the C++ version of the comdist function in the picante package
#'re22=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude_conspecifics = TRUE)
#'#test the result
#'expect_equal(as.matrix(re21),as.matrix(re22))
#'re23=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = TRUE)
#'#test on the C++ version of the comdist function in the picante package
#'re24=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude_conspecifics = TRUE)
#'expect_equal(as.matrix(re13),as.matrix(re14))
#'
#'#test on empty community
#'sample1=phylocom$sample
#'sample1[1,]=0
#'sample1[3,]=0
#'re=comdist(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE)
#'re2=comdist_C(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE)
#'expect_equal(re,re2)
#'re=comdist(sample1,cophenetic(phylocom$phylo),abundance.weighted=FALSE)
#'re2=comdist_C(sample1,cophenetic(phylocom$phylo),abundance.weighted=FALSE)
#'expect_equal(re,re2)
#'
#'re=comdistnt(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE)
#'re2=comdistnt_C(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE)
#'expect_equal(as.matrix(re),as.matrix(re2))
#'
#'#test the selection of community pairs
#'sample1=phylocom$sample
#'cal_pairs=matrix(TRUE,dim(sample1)[1],dim(sample1)[1])
#'cal_pairs[1,2]=FALSE
#'cal_pairs[2,1]=FALSE
#'re2=comdist_C(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE,cal_pairs)
#'expect_true( all(is.na(as.matrix(re2)[!cal_pairs])))
#'
#'re3=comdistnt_C(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE,TRUE,cal_pairs)
#'expect_true( all(is.na(as.matrix(re3)[!cal_pairs])))
#'
#'
#'#test the speed advantage
#'com=rCom(12000,10,win=owin(c(0,100),c(0,100)),ab="physignal",phy=list(br=rexp,phylosignal=1000))
#'com=quadratize(com,20,20)
#'sample=table(com$traits$ploti,com$traits$species)
#'dist=cophenetic(com$phylo)
#'t1=system.time(comdist(sample, dist, abundance.weighted=TRUE))
#'t2=system.time(comdist_C(sample, dist, abundance.weighted=TRUE))
#'expect_true(as.logical(t1[1]>t2[1]*100))
#'
#'t3=system.time(comdistnt(sample, dist, abundance.weighted=TRUE))
#'t4=system.time(comdistnt_C(sample, dist, abundance.weighted=TRUE))
#'expect_true(as.logical(t3[1]>t4[1]*100))
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