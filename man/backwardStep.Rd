\name{backwardStep}
\alias{backwardStep}
\title{To automatically select a best fitted clustering model based on pvalue}
\usage{
  backwardStep(fittedmodel, siglevel = 0.05)
}
\arguments{
  \item{fittedmodel}{a fm object represent a fitted
  clustering model for a point pattern}

  \item{siglevel}{a marginal significance level used to
  decide which part of model is significance.}
}
\description{
  To automatically select a best fitted clustering model
  based on pvalue
}
\details{
  The \code{backwardStep} accept a \code{fm} object from
  \code{\link{fitCluster}}, then check whether the trend of
  the model can be improved. Specifically, the improvement
  of trend is carried out by deleting insignificant
  environmental variables until all of then are significant
  or none of them left. significant pvalues of each
  environmental variables were estimated by considering any
  possible significant clustering in the residual.
}
\seealso{
  \code{\link{fitCluster}}, \code{\link{updateCluster}},
  \code{\link{varDecomp}}
}

