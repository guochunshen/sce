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
re1=comdist(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
#test on the C++ version of the comdist function in the picante package
re2=comdist_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
#test the result
expect_equal(re1,re2)

#test the speed advantage
com=rCom(12000,10,win=owin(c(0,100),c(0,100)),ab="physignal",phy=list(br=rexp,phylosignal=1000))
com=quadratize(com,20,20)
sample=table(com$traits$ploti,com$traits$species)
dist=cophenetic(com$phylo)
t1=system.time(comdist(sample, dist, abundance.weighted=TRUE))
t2=system.time(comdist_C(sample, dist, abundance.weighted=TRUE))
expect_true(as.logical(t1[1]>t2[1]*100))
}

