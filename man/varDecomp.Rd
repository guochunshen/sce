\name{varDecomp}
\alias{varDecomp}
\title{To estimate the proportion of variance explained by habitat, internal cluster and random process from a given fitted cluster model}
\usage{
  varDecomp(fittedmodel, r = c(0:80), R = 10, delta = 1,
    simple = FALSE)
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

  \item{simple}{if it is true, only PVH and variance of
  habitat and spatial autocorrelation were calculated}
}
\value{
  a vector of length 6 contains
  \code{"PVH","PVR","varHabitat","varInterCluster","varPoi","varNonrand"}.
  Specifically, "PVH": proportion of variance explained by
  Habitat, thus proportion of variance explained by
  internal clustering processes is \bold{PVI}=1-PVH "PVR":
  proportion of variance explained by poisson Random
  process, thus proportion of variance explained by the
  nonrandom processes (\bold{PVNR}) is 1-PVR "varHabitat"
  and "varInterCluster" are the variance explained by
  habtiat and internal cluster process in the log intensity
  function of the point pattern; "varPoi" and "varNonrand"
  are the variances explained by poisson random and
  nonrandom processes in the intensity function of the
  point pattern.
}
\description{
  To estimate the proportion of variance explained by
  habitat, internal cluster and random process from a given
  fitted cluster model
}

