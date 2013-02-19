\name{subset.scp}
\alias{subset.scp}
\title{Extract a subject set of a community}
\usage{
  subset(x,i) \cr x[i]
}
\arguments{
  \item{x}{an scp object}

  \item{i}{in numeric or logical index of needed data}
}
\description{
  Extract a subject set of a community
}
\examples{
data(testData)

#get the first one hundred individuals
testData[1:100]
#or
subset(testData,1:100)

#get a perticular species
sel=testData$traits$species==testData$sp[2]
testData[sel]
#or
subset(testData,sel)
}

