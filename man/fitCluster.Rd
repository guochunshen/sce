\name{fitCluster}
\alias{fitCluster}
\title{Fitting the each group individual mapped population in the community by a clustering model}
\usage{
  fitCluster(com, trend = ~1, cluster = "LGCP",
    sigTest = FALSE,
    ctlpars = list(rmax = 25, rmin = 3, bw = 5, sigma2 = 3, alpha = 10, nurange = c(Inf, 0.5), edgecor = "translate", nsim = 10, r = seq(0, 60, 2), siglevel = 0.05),
    ...)
}
\arguments{
  \item{com}{a \link{scp} object of a given community}

  \item{trend}{a specific formular for modeling trend.}

  \item{cluster}{define how internal cluster assembled with
  intensity function. possible choose are "poisson" and
  "LGCP"}

  \item{sigTest}{a logical flag to test significance of
  extra aggregation in residual and habitat variables}

  \item{ctlpars}{control parameters used in model fitting.
  See meaning of each par in Details.}
}
\value{
  a list of fm object contains fitted parameters for each
  group
}
\description{
  \code{fitCluster} fit the point pattern of each groups
  (species, DBH class) by a homogeneous/heterogeneous
  Cox/cluster process.
}
\details{
  Three steps were carried out: 1. fit the point pattern by
  a heterogeneous poisson process; 2. calculate pair
  correlation function of the residual (using adaptive
  estimation, see details in \code{\link{pcf_adaptive}}).
  3. estimate internal clustering parameters by minimum
  contrast method.

  the \code{ctlpars} contains parameter values that are
  very sensitive to the final result. it controls
  calculations of the pari correlation function, the
  minimum contrast method and significance tests of habitat
  and extra internal aggregation. Specifically, \code{rmax}
  controls the spatial range that calculate the pair
  correlation function in the \code{\link{pcf_adaptive}}
  function. \code{rmax} and \code{rmin} also define the
  range of distance considered in the minimum contrast
  method n estimating parameters of internal aggregation.
  \code{nurange} controls the sharp of marten covariance
  function in the minimum control method. \code{bw}
  controls the bindwidth of kernel estimation in
  \code{\link{pcf_adaptive}}. \code{sigma2}, \code{alpha}
  and \code{edgecor} are initial parameters of the fitted
  clustering process.  \code{nsim}, \code{r} and
  \code{edgecor} are parameters used in testing
  significance of aggregated pattern in the residual. it
  adapted a goodness-of-fit test method. Thus they control
  the number of simulations, interpoint ditance range and
  edge correction method respectively.
}
\examples{
#load the testData set
data(testData)

#remove rare species
com=removeRareSpecies(testData,600)

#fit pattern of each species by a best cluster model
fittedmodel=fitCluster(com,~elev+grad,sigTest=TRUE)
}
\seealso{
  \link{sigHabitatTest},
  \link{sigAggreResidualTest},\link{updateCluster}
  ,\link{varDecomp}
}

