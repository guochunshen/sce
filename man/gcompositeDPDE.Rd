\name{gcompositeDPDE}
\alias{gcompositeDPDE}
\title{The general composite distance based population density estimator (DPDE)}
\usage{
  gcompositeDPDE(ptoe_r, etoe_r, k, area)
}
\arguments{
  \item{ptoe_r}{neighborhood distances matrix from the
  focal point/location to its nearest neighbor, the number
  of column represent the number of equal angle sector used
  in sampling distances}

  \item{etoe_r}{neighborhood distances matrix from the
  focal event to its nearest neighbor}

  \item{k}{the order of neighbor}

  \item{area}{the area in which population size want to be
  estimated}
}
\description{
  The general composite distance based population density
  estimator (DPDE)
}
\examples{
ptoe_r=matrix(runif(100,0,4),nrow=50,ncol=2)
#ptoe_r and etoe_r can have different sample sizes
etoe_r=matrix(runif(120,0,4),nrow=60,ncol=2)

gcompositeDPDE(ptoe_r,etoe_r,k=2,area=100)
}

