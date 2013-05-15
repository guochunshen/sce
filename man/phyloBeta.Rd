\name{phyloBeta}
\alias{phyloBeta}
\title{calculate the the quadrat based phylogenetic beta diversity with significant test by simulations}
\usage{
  phyloBeta(com, phyd, Fun, nsim, rmax = 50, alpha = 0.05,
    ...)
}
\arguments{
  \item{com}{a quadratized scp object trait for each
  individual and phylogeney}

  \item{phyd}{phylogenetic distance matrix}

  \item{nsim}{number of shuffling tip from the given
  phylongeny}

  \item{rmax}{maximum spatial distance between two
  communities}

  \item{...}{other parameters past to
  \code{\link{comdist_C}} function.}
}
\description{
  calculate the the quadrat based phylogenetic beta
  diversity with significant test by simulations
}
\examples{
library(testthat)

com=rCom(1000,10,win=owin(c(0,100),c(0,100)),ab="physignal",phy=list(br=rexp,phylosignal=1000))
com=quadratize(com,20,20)
phyd=cophenetic(com$phylo)
re1=phyloBeta(com,phyd=phyd,Fun=comdist_C,nsim=10,abundance.weighted=FALSE)
re2=phyloBeta(com,phyd=phyd,Fun=comdist_C,nsim=10,abundance.weighted=TRUE)

expect_true(any(re1$real!=re2$real))

re1=phyloBeta(com,phyd=phyd,Fun=comdistnt_C,nsim=10,abundance.weighted=FALSE)
re2=phyloBeta(com,phyd=phyd,Fun=comdistnt_C,nsim=10,abundance.weighted=TRUE)

expect_true(any(re1$real!=re2$real))

expect_true(max(re1$r)<50)
}

