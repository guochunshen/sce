\name{comdistnt_C}
\alias{comdistnt_C}
\title{An C++ implementation of inter-community mean nearest taxon distance}
\usage{
  comdistnt_C(comm, dis, abundance.weighted = FALSE,
    exclude.conspecifics = FALSE, cal_pairs = NULL)
}
\arguments{
  \item{comm}{Community data matrix}

  \item{interspecific}{distance matrix}

  \item{abundance.weighted}{Should mean nearest taxon
  distances from each species to species in the other
  community be weighted by species abundance? (default =
  FALSE)}

  \item{exclude.conspecifics}{Should conspecific taxa in
  different communities be exclude from MNTD calculations?
  (default = FALSE)}
}
\description{
  An C++ implementation of inter-community mean nearest
  taxon distance
}
\note{
  it is an C++ implementation of the comdistnt function in
  picante package
}
\examples{
library(testthat)

re21=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = FALSE)
#test on the C++ version of the comdist function in the picante package
re22=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = FALSE)
#test the result
expect_equal(as.matrix(re21),as.matrix(re22))
re23=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = FALSE)
#test on the C++ version of the comdist function in the picante package
re24=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = FALSE)
expect_equal(as.matrix(re13),as.matrix(re14))

re21=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = TRUE)
#test on the C++ version of the comdist function in the picante package
re22=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE, exclude.conspecifics = TRUE)
#test the result
expect_equal(as.matrix(re21),as.matrix(re22))
re23=comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = TRUE)
#test on the C++ version of the comdist function in the picante package
re24=comdistnt_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE, exclude.conspecifics = TRUE)
expect_equal(as.matrix(re13),as.matrix(re14))
}

