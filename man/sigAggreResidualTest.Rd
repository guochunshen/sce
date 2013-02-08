\name{sigAggreResidualTest}
\alias{sigAggreResidualTest}
\title{To test whether there is any aggregative residual existed in the unexplianed pattern by habitat}
\usage{
  sigAggreResidualTest(fittedmodel, nsim = 10,
    r = seq(0, 60, 2), edcor = "translate")
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

