\name{updateCluster}
\alias{updateCluster}
\title{To update a model manually or automatically}
\usage{
  updateCluster(fittedmodel, trend = NULL, siglevel = 0.05)
}
\arguments{
  \item{fittedmodel}{a fm object represent a fitted
  clustering model for a point pattern}

  \item{trend}{a \code{\link{formula}} object represent the
  structure of the intensity trend. if it is null, the most
  insignificant habitat will be removed from the formula.}

  \item{siglevel}{a marginal significance level used to
  decide which part of model is significance.}
}
\description{
  To update a model manually or automatically
}

