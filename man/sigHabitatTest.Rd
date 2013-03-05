\name{sigHabitatTest}
\alias{sigHabitatTest}
\title{To calculate pvalues of habitat varaibles under consideration of extra aggregation existed in residual.}
\usage{
  sigHabitatTest(fittedModel, clusterResidual = TRUE)
}
\arguments{
  \item{fittedModel}{a \link{fm} object representing a
  fitted cluster point process for the given point
  pattern.}

  \item{clusterResidual}{a logical flag to consider any
  left cluster residual if it is true.}
}
\description{
  To calculate pvalues of habitat varaibles under
  consideration of extra aggregation existed in residual.
}
\details{
  if there is significant clustered residual
  (\code{clusterResidual=TRUE}), significance of the
  habitat will be evaludated by considering those extra
  aggregation in residual. otherwise
  (\code{clusterResidual=FALSE}), regular significant test
  method of regression coefficients was carried out.
}
\examples{
data(testData)

sp1=subset(testData,testData$traits$species=="ACALDI")

fm=fitCluster(sp1,~elev+grad,sigTest=FALSE)

sigHabitatTest(fm)
}

