\name{rPopNeg}
\alias{rPopNeg}
\title{Generate simulated population with negative binomial abundance distribution}
\usage{
  rPopNeg(lam, a, win, nxy)
}
\arguments{
  \item{lam}{expected individual intensity of a population}

  \item{a}{the aggregated parameter in the negative
  binomial distribution. the smaller the a is, the more
  aggregate the populaiton is.}

  \item{win}{a owin object represent the range of the
  population}

  \item{nxy}{a vector of number of quadrats in the x and y
  axis respectively.}
}
\value{
  a scp object
}
\description{
  Generate simulated population with negative binomial
  abundance distribution
}

