\name{rCommunity}
\alias{rCommunity}
\title{create a community driven by single process like dispersal limitation, habitat association or competition}
\usage{
  rCommunity(com_control)
}
\arguments{
  \item{com_control}{a list object contains expected
  properties about the simulated community.}
}
\value{
  a list contains: a marked point pattern and a ecological
  difference matrix
}
\description{
  create a community driven by single process like
  dispersal limitation, habitat association or competition
}
\details{
  The com_control has follow field about the community:

  N: total number individuals in the community

  S: total number species in the community

  win: an \code{\link{owin}} object represent the spatial
  range of the community

  ab: a character discrebing the species abundance
  distribution of the target community. currently avaible
  options are "uniform", "logseries" and "lognormal"

  intra: a character representing intraspecific
  relationship. currently avaible options are "Poisson",
  "dispersal"

  inter: a string representing interspecific relationship.
  currently avaible options are "independent" and
  "dependent"

  covr: a logical flag to add species habitat association
  into the community. currently, habitat is simpled modeled
  by the sine function.  it always be FALSE in community
  assembled by competition.

  covrscale: the spatial scale of the covariance.

  beta: a real number defines the competition intensity
  between species. the large the number is, the intense the
  competition is.

  r: a real distance value. competition only happend within
  this distance.

  R2: define the strength of correlation between species
  character and their phylogenetic relationship
}
\examples{
com_control=list(S=10,N=1000,ab="lognormal",win=owin(c(0,100),c(0,100)),
intra="Poisson",inter="independent",R2=0.9,covr=FALSE,beta=2,r=5,covrscale=16)

#Poisson distribution
rCommunity(com_control)

#dispersal
com_control$intra="dispersal"
rCommunity(com_control)

#habitat
com_control$intra="Poisson"
com_control$covr=TRUE
rCommunity(com_control)

#compeition, relatively time consuming
com_control$covr=FALSE
com_control$inter="dependent"
rCommunity(com_control)
}

