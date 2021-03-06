\name{rmypoispp}
\alias{rmypoispp}
\title{generate a poisson distribution point pattern, avoid 0 and very large number of point}
\usage{
  rmypoispp(lambda, Nexpect)
}
\arguments{
  \item{lambda}{an im object represent the intensity map of
  the expected point pattern}

  \item{Nexpect}{the expect number of points in the
  observed window}
}
\description{
  generate a poisson distribution point pattern, avoid 0
  and very large number of point
}
\details{
  This is a warp function of \code{\link{rpoispp}} in
  spatstat package, and add two limitations into it. one is
  that avoid generating too much data points that commonly
  generated by a wrong lambda. thus it avoid the memory
  problem. The second improve is elimate the possibility of
  generating a point pattern with 0 point. this commonly
  happened in the case of very small lambda. although there
  is no problem with the zero point pattern, it will cause
  any other unexpected problem in other functions that use
  the simulated data point.
}

