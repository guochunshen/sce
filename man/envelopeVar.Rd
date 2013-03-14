\name{envelopeVar}
\alias{envelopeVar}
\title{calculate envelope of variances explained different fractions of the model}
\usage{
  envelopeVar(fittedmodel, nsim = 9, conf_level = 0.95,
    r = c(0:80), R = 10, delta = 1, simple = FALSE)
}
\arguments{
  \item{fittedmodel}{a fm object returned by the
  \code{\link{fittedmodel}} function}

  \item{nsim}{number of simulations used to generate the
  envelope}

  \item{conf_level}{significance level of the confidence
  interval}

  \item{r,R,delta}{parameters used in the
  \code{\link{varDecomp}} function}

  \item{simple}{a logical flag in \code{\link{varDecomp}}
  to calculate the PVH only if it is true}
}
\description{
  calculate envelope of variances explained different
  fractions of the model
}

