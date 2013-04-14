\name{fastBM}
\alias{fastBM}
\title{Fast Brownian simulation}
\usage{
  fastBM(tree, a = 0, mu = 0, sig2 = 1,
    bounds = c(-Inf, Inf), internal = FALSE, nsim = 1)
}
\arguments{
  \item{tree}{is a phylogenetic tree in "phylo" format.}

  \item{a}{a value for ancestral state at the root node.}

  \item{mu}{an optional value for the mean of random normal
  changes along branches of the tree - can be used to
  simulate a trend if mu!=0.}

  \item{sig2}{instantaneous variance of the BM process.}

  \item{bounds}{a vector with the lower and upper bounds
  (respectively) for bounded Brownian simulation - by
  default simulation is unbounded.}

  \item{internal}{logical value indicating whether or not
  to return states for internal nodes.}

  \item{nsim}{number of simulations.}
}
\description{
  Fast Brownian simulation
}
\details{
  it is an only copy of fastBM from phytools package.
  please see details in that package.
}

