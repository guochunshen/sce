\name{phyMarkCorr}
\alias{phyMarkCorr}
\title{The phylogenetic mark correlation function with species shuffling null model}
\usage{
  phyMarkCorr(com, phyd, focalij = NULL, rmax, step = 1,
    nsim = 19, alpha = 0.05, scale = TRUE)
}
\arguments{
  \item{com}{a scp object}

  \item{phyd}{a phylogenetic distance matrix with species
  colnames and colnames}

  \item{focalij}{a logical vector indicate which
  individuals are focal individuals.  if it is null, all
  individuals are the focal points.}

  \item{rmax}{the maximum distance to calculate the
  phylogenetic mark correlatin function}

  \item{step}{the step of distance from 0 to quantify
  spaital phylogenetic structure}

  \item{nsim}{number of shuffling tip labels in the species
  shuffling null model}

  \item{alpha}{sinificant level}

  \item{scale}{a logical to indicate wether the
  phylogenetic mark correlation function should be scaled.}
}
\description{
  The phylogenetic mark correlation function with species
  shuffling null model
}
\details{
  explained detailed further
}
\examples{
data(testData)

phyMarkCorr(testData,testData$phylo,rmax=20)
}

