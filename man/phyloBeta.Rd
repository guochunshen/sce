\name{phyloBeta}
\alias{phyloBeta}
\title{calculate the the quadrat based phylogenetic beta diversity with significant test by simulations}
\usage{
  phyloBeta(com, Fun, nsim, rmax = 50, alpha = 0.05, ...)
}
\arguments{
  \item{com}{a quadratized scp object trait for each
  individual and phylogeney}

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

