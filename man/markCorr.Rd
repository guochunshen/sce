\name{markCorr}
\alias{markCorr}
\title{The c++ implemented mark correlation function with mark random labeling null model}
\usage{
  markCorr(com, markName, r, testFunName, nsim = 10,
    h = 0.5, exclude_conspecific = FALSE, normalize = TRUE,
    alpha = 0.05)
}
\arguments{
  \item{com}{a scp object}

  \item{markName}{the name of the mark contained in
  \code{com} to calculate the mark pattern}

  \item{r}{a vector of neighborhood distances to
  quantifying the mark correlation function}

  \item{testFunName}{name of the test function, see details
  of the current implemented test functions.}

  \item{nsim}{number of shuffling the values of mark}

  \item{h}{the torrent bin.}

  \item{exclude_conspecific}{a logical flag to whether
  include the conspecific indiviual or not}

  \item{normalize}{a logical flag to normalize the mark
  correlation function}

  \item{alpha}{significant level used to calculate
  pointwise confidence interval}
}
\description{
  The c++ implemented mark correlation function with mark
  random labeling null model
}
\details{
  it calculate the mark correlation function for individual
  based numerical marks, e.g. DBH. current implemented test
  functions are

  testFunName: abdif (absolute difference): abs(m1-m2)
  testFunName: diff (difference from focal to neighbor):
  m1-m2 testFunName: sum (sum of the marks of two points):
  m1+m2 testFunName: reldif (relative absolute difference):
  abs(m1-m2)/(m1+m2)
}
\examples{
data(testData)
re1=markCorr(com=testData,markName="dbh",r=1:10,testFunName="tf1",nsim=10,h=1,exclude_conspecific=TRUE)
}

