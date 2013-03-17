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
  pari correlation function.
}
\examples{
data(testData)

sp1=subset(testData,testData$traits$species=="ACALDI")

fm=fitCluster(sp1,~elev+grad)

sigAggreResidualTest(fm,nsim=10,r=0:30)
}

