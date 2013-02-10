\name{varDecomp}
\alias{varDecomp}
\title{To estimate the proportion of variance explained by habitat, internal cluster and random process from a given fitted cluster model}
\usage{
  varDecomp(fittedmodel, r = c(0:80), R = 10, delta = 1)
}
\arguments{
  \item{fittedmodel}{a \code{\link{fm}} object representing
  a fitted clustering model}

  \item{r}{spatial range used in pair correlation function
  to caluclate the proportaion of variance}

  \item{R}{the quadrat scale used to evaluate the variance
  of poisson noise}

  \item{delta}{use delta x delta grid for numerical
  integration.}
}
\description{
  To estimate the proportion of variance explained by
  habitat, internal cluster and random process from a given
  fitted cluster model
}

