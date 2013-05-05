\name{gsimpleDPDE}
\alias{gsimpleDPDE}
\title{The general simple distance based population density estimator (DPDE)}
\usage{
  gsimpleDPDE(r_samples, k, area)
}
\arguments{
  \item{r_samples}{a matrix of distance, distances in each
  column comes from one section of q equal-angle sectors.}

  \item{k}{the order of nearest neighbor}

  \item{area}{the area of region that want to estimate
  number of individuals}
}
\description{
  The general simple distance based population density
  estimator (DPDE)
}
\examples{
r_samples=matrix(runif(100,0,4),nrow=50,ncol=2)
gsimpleDPDE(r_samples,k=1,area=100)
}

