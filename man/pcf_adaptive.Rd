\name{pcf_adaptive}
\alias{pcf_adaptive}
\title{estimate the inhomogeneous pair correlation function by adaptive method}
\usage{
  pcf_adaptive(data.ppp, lambda, maxt, bw, adaptive = T,
    kerneltype = 1)
}
\arguments{
  \item{data}{a \code{spatstat ppp} point patttern object}

  \item{lambda}{vector containing intensity for each point
  in data}

  \item{maxt}{maximal distance t for which g(t) is
  estimated}

  \item{bw}{band with for kernel estimation}

  \item{adaptive}{if true, Snethlage (2001)'s adaptive
  estimate is used}

  \item{kerneltype}{if 1, the uniform kernel is used,}
}
\value{
  an \code{fv} object
}
\description{
  \code{pcf_adaptive} is only a warp function of
  \link{mygest} function in InhomCluster package to fit the
  formate of fv object defined in spatstat package.
}
\details{
  This function computes a kernel estimate of the
  g-function using the translation invariant edge
  correction (Ohser, 1983). The procedure only works for
  point patterns observed within a rectangular observation
  window.
}

