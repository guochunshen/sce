\name{varDecomp}
\alias{varDecomp}
\title{fit the given point pattern by a best cluster model
and then estimate the proportion of variance explained by habitat and internal cluster}
\usage{
  varDecomp(com, multicore = FALSE, mc.cores = 1,
    pars = list(rmax = 70, rmin = 5, bw = 15, nsim = 39, siglevel = 0.05, sigma2 = 3, alpha = 10, nurange = c(Inf, 0.5), q = 2, siglevel2 = 0.05, agtest = TRUE, edgecor = "translate"))
}
\arguments{
  \item{com}{a \link{scp} object that contains distribution
  of trees and environmental variables}

  \item{multicore}{a logical flag indicate whether using
  multicore to speedup calculation.  physical support of
  multicore computation is needed first}

  \item{mc.cores}{number of cores used. it also means how
  many species were under calculation simutatnously.}

  \item{pars}{sensitive internal parameters related to
  model fitting and significance test.  do not change it
  unless what you are doing.}
}
\description{
  \code{varDecomp} using three steps to estimate proportion
  of variance in spatial point pattern. First, it fit the
  given point pattern by a cluster model, in which
  significance of each habitat variables, signficance of
  remaining cluster were estimated. Then a best cluster
  model was chosen by the backward stepwise selection
  method.
}

