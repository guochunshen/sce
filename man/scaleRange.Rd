\name{scaleRange}
\alias{scaleRange}
\title{scale a vector's range and keep the relationship between elements}
\usage{
  scaleRange(x, min, max)
}
\arguments{
  \item{x}{a numeric vector or maxtrix to be scale the
  range}

  \item{min}{the minimum value of the output vector}

  \item{max}{the maxmum value of the output vector}
}
\description{
  scale a vector's range and keep the relationship between
  elements
}
\note{
  every elements in x should be positive
}
\examples{
x=runif(100,1,100)
x2=scaleRange(x,0.5,1)
cor(x,x2)

x=matrix(x,10,10)
x3=scaleRange(x,0.5,1)
cor(as.numeric(x),as.numeric(x3))
}

