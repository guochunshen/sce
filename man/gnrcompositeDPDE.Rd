\name{gnrcompositeDPDE}
\alias{gnrcompositeDPDE}
\title{The general composite DPDE based on negative binomial distribution}
\usage{
  gnrcompositeDPDE(ptoe_r, etoe_r, k, area)
}
\description{
  The general composite DPDE based on negative binomial
  distribution
}
\examples{
ptoe_r=matrix(runif(100,1,4),nrow=50,ncol=2)
etoe_r=matrix(runif(100,1,3),nrow=50,ncol=2)

gnrcompositeDPDE(ptoe_r,etoe_r,k=2,area=100)
}

