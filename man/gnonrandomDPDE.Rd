\name{gnonrandomDPDE}
\alias{gnonrandomDPDE}
\title{The general nonrandom distance based population density estimator (DPDE)}
\usage{
  gnonrandomDPDE(r_samples, dtype = "ptoe", k, area)
}
\arguments{
  \item{r_samples}{neighborhood distance matrix. The number
  of column represents the number of equal angle sectors
  used in distance sampling}

  \item{dtype}{type of neighborhood distance in
  \code{r_sample}. it could be "ptoe" (point to event
  distance) or "etoe" (event to event distance).}

  \item{k}{the order of neighborhood distance}

  \item{area}{area the area in which population size want
  to be estimated}
}
\description{
  The general nonrandom distance based population density
  estimator (DPDE)
}
\details{
  This population density esimator is based on a assumption
  of population density distribution across quadrat
  following a negative binomial distribution. In this case,
  the forms of estimator based on point to event distances
  and event to event distances are different.
}
\examples{
r_samples=matrix(runif(100,0,4),nrow=50,ncol=2)
gnonrandomDPDE(r_samples,dtype="ptoe",k=1,area=100)
gnonrandomDPDE(r_samples,dtype="etoe",k=1,area=100)
}

