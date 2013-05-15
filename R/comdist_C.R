#'
#' A c++ implementation of the inter-community mean pairwise distance
#' 
#'@param comm Community data matrix, rows are the sample plot, cols are the species name
#'@param dis interspecific distance matrix
#'@param abundance.weighted Should mean pairwise distances separating species in two communities be weighted by species abundances? (default = FALSE)
#'@param cal_pairs a logical matrix to indicate which pair of communities should be used to calculate the index
#'
#'@note this is an c++ implementation of comdist function in the picanate package.
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


comdist_C<-function (comm, dis, abundance.weighted = FALSE, cal_pairs=NULL) 
{
  x <- as.matrix(comm)
  dat <- match.comm.dist(comm, dis)
  x <- dat$comm
  dis <- as.matrix(dat$dist)
  if (!abundance.weighted) {
    x <- decostand(x, method = "pa")
  }
  N <- dim(x)[1]
  S <- dim(x)[2]
  x <- decostand(x, method = "total", MARGIN = 1)
  if(is.null(cal_pairs)){
    cal_pairs=matrix(TRUE,nrow=N,ncol=N)
  }
  
  comdist=rep(0,N^2)
  comdist=.C("comdist",as.integer(N),as.integer(S),as.double(dis),as.double(x),as.integer(cal_pairs),as.double(comdist))[[6]]
  dim(comdist)=c(N,N)
  #set the uninteresting community pair equals to NA
  
  comdist[!cal_pairs]=NA
  row.names(comdist) <- row.names(x)
  colnames(comdist) <- row.names(x)
  return(as.dist(comdist))
}