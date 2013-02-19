\name{subsetScp}
\alias{subsetScp}
\title{Extract a subject set of a community}
\usage{
  subsetScp(x, i)
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
subsetScp(testData,1:100)

#get a perticular species
sel=testData$traits$species==testData$sp[2]
subsetScp(testData,sel)
}

