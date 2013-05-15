\name{comdist_C}
\alias{comdist_C}
\title{A c++ implementation of the inter-community mean pairwise distance}
\usage{
  comdist_C(comm, dis, abundance.weighted = FALSE,
    cal_pairs = NULL)
}
\arguments{
  \item{comm}{Community data matrix, rows are the sample
  plot, cols are the species name}

  \item{dis}{interspecific distance matrix}

  \item{abundance.weighted}{Should mean pairwise distances
  separating species in two communities be weighted by
  species abundances? (default = FALSE)}

  \item{cal_pairs}{a logical matrix to indicate which pair
  of communities should be used to calculate the index}
}
\description{
  A c++ implementation of the inter-community mean pairwise
  distance
}
\note{
  this is an c++ implementation of comdist function in the
  picanate package.
}
\examples{
library(testthat)

data(phylocom)
re11=comdist(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
#test on the C++ version of the comdist function in the picante package
re12=comdist_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
#test the result
expect_equal(re11,re12)
re13=comdist(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE)
#test on the C++ version of the comdist function in the picante package
re14=comdist_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE)
expect_equal(re13,re14)

re21=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = FALSE)
#test on the C++ version of the comdist function in the picante package
re22=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude_conspecifics = FALSE)
#test the result
expect_equal(as.matrix(re21),as.matrix(re22))
re23=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = FALSE)
#test on the C++ version of the comdist function in the picante package
re24=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude_conspecifics = FALSE)
expect_equal(as.matrix(re13),as.matrix(re14))

re21=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = TRUE)
#test on the C++ version of the comdist function in the picante package
re22=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude_conspecifics = TRUE)
#test the result
expect_equal(as.matrix(re21),as.matrix(re22))
re23=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = TRUE)
#test on the C++ version of the comdist function in the picante package
re24=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude_conspecifics = TRUE)
expect_equal(as.matrix(re13),as.matrix(re14))

#test on empty community
sample1=phylocom$sample
sample1[1,]=0
sample1[3,]=0
re=comdist(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE)
re2=comdist_C(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE)
expect_equal(re,re2)
re=comdist(sample1,cophenetic(phylocom$phylo),abundance.weighted=FALSE)
re2=comdist_C(sample1,cophenetic(phylocom$phylo),abundance.weighted=FALSE)
expect_equal(re,re2)

re=comdistnt(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE)
re2=comdistnt_C(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE)
expect_equal(as.matrix(re),as.matrix(re2))

#test the selection of community pairs
sample1=phylocom$sample
cal_pairs=matrix(TRUE,dim(sample1)[1],dim(sample1)[1])
cal_pairs[1,2]=FALSE
cal_pairs[2,1]=FALSE
re2=comdist_C(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE,cal_pairs)
expect_true( all(is.na(as.matrix(re2)[!cal_pairs])))

re3=comdistnt_C(sample1,cophenetic(phylocom$phylo),abundance.weighted=TRUE,TRUE,cal_pairs)
expect_true( all(is.na(as.matrix(re3)[!cal_pairs])))


#test the speed advantage
com=rCom(12000,10,win=owin(c(0,100),c(0,100)),ab="physignal",phy=list(br=rexp,phylosignal=1000))
com=quadratize(com,20,20)
sample=table(com$traits$ploti,com$traits$species)
dist=cophenetic(com$phylo)
t1=system.time(comdist(sample, dist, abundance.weighted=TRUE))
t2=system.time(comdist_C(sample, dist, abundance.weighted=TRUE))
expect_true(as.logical(t1[1]>t2[1]*100))

t3=system.time(comdistnt(sample, dist, abundance.weighted=TRUE))
t4=system.time(comdistnt_C(sample, dist, abundance.weighted=TRUE))
expect_true(as.logical(t3[1]>t4[1]*100))
}

