\name{gofTest}
\alias{gofTest}
\title{a general goodness-of-fit (gof) test on the fitted cluster model}
\usage{
  gofTest(fittedmodel, SFun, rRange = c(0, 10), nsim = 10,
    r = seq(0, 20, 1), ...)
}
\arguments{
  \item{a}{fittedmodel returned by \code{\link{fitCluster}}
  function}

  \item{SFun}{a summary statistic of spaital point pattern}

  \item{rRange}{the neighborhood distance range used to
  check performance of the model.}

  \item{nsim}{number of simulated patterns.}

  \item{r}{Numeric vector. The values of the argument r at
  which summary statistic should be evaluated.}

  \item{...}{other parameters passed to the summary
  statistic}
}
\value{
  a vector of pvalue, each pvalue corresponding to a edge
  correction method.
}
\description{
  a general goodness-of-fit (gof) test on the fitted
  cluster model
}
\details{
  We could use a similar technique to do a goodness-of-fit
  test for the fitted log Gaussian Cox process. Consider
  e.g. the F-function. We can compute an estimate F_dat
  from the data and we can obtain estimated F-functions
  F_i, i=1,...,n from simulations of the fitted log
  Gaussian Cox process. We can then approximate the
  theoretical value of the F function by the mean F_mean of
  the simulated F_i, i=1,...,n. Our goodness of fit
  statistic then becomes

  D=integral from r_min to r_max of the squared distance
  (F_dat(r)-F_mean(r))^2.

  Similarly we obtain simulated values of D as

  D_i=integral from r_min to r_max of the squared distance
  (F_i(r)-F_mean(r))^2.

  Finally our goodness of fit p-value is (m+1)/(n+1) where
  m is the number of D_i which are larger than D.

  The above procedure can be used with all the summary
  statistics. However the power of the goodness of fit test
  is weak with g since g was used for the model fitting
  (hence we should expect the theoretical g to be close to
  the g estimated from the data). Since K is just the
  integrated g, the goodness of fit test will also be weak
  with K. Nevertheless, we could carry out the goodness of
  fit test both with g, K, F, G, and J.

  All the statistics are available in spatstat. Note
  however, that "inhomogeneous" versions of F, K, and G are
  not available so here we have to use the "stationary"
  versions.

  an insignificant pvalue means the fitted model has no
  significant different with the observed data.
}
\examples{
#load one species data
data(testData)
sp1=subset(testData,testData$traits$species=="BEILPE")

#fit a cluster model
fittedmodel=fitCluster(sp1,~elev+grad)

#goodness-of-fit test
gofTest(fittedmodel,SFun=Fest,rRange=c(0,10),nsim=20,r=seq(0,20,1))
}

