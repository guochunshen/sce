\name{quadratize}
\alias{quadratize}
\title{divide a large plot into small quadrats}
\usage{
  quadratize(com, nx, ny)
}
\arguments{
  \item{com}{a scp object}

  \item{nx}{number of quadrat in the x axis}

  \item{ny}{number of qudart in the y axis}
}
\value{
  a new scp object with a "ploti" trait in the traits of
  the community it also has an attribute "spaced" contains
  spatial distance between quadrats
}
\description{
  divide a large plot into small quadrats
}
\examples{
data(testData)

newdata=quadratize(testData,5,5)

table(newdata$traits$ploti)
}

