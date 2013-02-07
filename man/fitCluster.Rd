\name{fitCluster}
\alias{fitCluster}
\title{Fitting the each group individual mapped population in the community by a clustering model}
\usage{
  fitCluster(com, trend = ~1, cluster = "LGCP",
    covmodel = list(model = "matern", nu = 0.3),
    group = NULL,
    ctlpars = list(rmax = 25, rmin = 3, bw = 5, sigma2 = 3, alpha = 10, nurange = c(Inf, 0.5), q = 2, edgecor = "translate"),
    ...)
}
\arguments{
  \item{com}{a \link{scp} object of a given community}

  \item{trend}{a specific formular for modeling trend.}

  \item{cluster}{define how internal cluster assembled with
  intensity function}

  \item{covmodel}{define the exact covariance function
  under LGCP model}

  \item{group}{factor, defination of unit of individuals
  used to fit model}

  \item{ctlpars}{control parameters used in model fitting.}
}
\value{
  a list of fitted parameters for each group
}
\description{
  \code{fitCluster} fit the point pattern of each groups
  (species, DBH class) by a homogeneous/heterogeneous
  Cox/cluster process.
}
\details{
  Three steps were carried out: 1. fit the point pattern by
  a heterogeneous poisson process; 2. calculate pair
  correlation function of the residual. 3. estimate
  internal clustering parameters by minimum contrast
  method.
}

