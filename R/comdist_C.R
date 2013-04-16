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
#'re1=comdist(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
#'#test on the C++ version of the comdist function in the picante package
#'re2=comdist_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
#'#test the result
#'expect_equal(re1,re2)
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
  comdist = commdistInner(N,dis,x,cal_pairs)
  #set the uninteresting community pair equals to NA
  comdist[!cal_pairs]=NA
  row.names(comdist) <- row.names(x)
  colnames(comdist) <- row.names(x)
  return(as.dist(comdist))
}