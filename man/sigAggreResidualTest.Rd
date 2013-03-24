\name{sigAggreResidualTest}
\alias{sigAggreResidualTest}
\title{To test whether there is any aggregative residual existed in the unexplianed pattern by habitat}
\usage{
  sigAggreResidualTest(fittedmodel, nsim = 10,
    r = seq(0, 60, 2), edcor = "translate", ntry = 10)
}
\arguments{
  \item{fittedmodel}{a \link{fm} object represented a
  fitted cluster point process model for a point pattern.}

  \item{nsim}{number of relazations used to simulate point
  pattern accroding to the model.}

  \item{r}{distance range used in goodness of fit test.}

  \item{edcor}{edge correction method. see details in
  \link{pcf}}
}
\description{
  To test whether there is any aggregative residual existed
  in the unexplianed pattern by habitat
}
\details{
  This is a goodness-of-fit test based on the heterogeneous
  pari correlation function. since
  \code{\link{pcf_adaptive}} is good at lower memory and
  worse on scater data point pattern (e.g. number of
  individuals below than 50). Thus we use
  \code{\link{pcf_adaptive}} on abundant species, and use
  \code{\link{pcf}} on rare species.
}
\examples{
data(testData)

sp1=subset(testData,testData$traits$species=="ACALDI")

fm=fitCluster(sp1,~elev+grad)

#significant residual
sigAggreResidualTest(fm,nsim=10,r=0:30)

sp2=subset(testData,testData$traits$species=="HIRTTR")

fm2=fitCluster(sp2,~elev+grad)

sigAggreResidualTest(fm2,nsim=10,r=0:30)
}

