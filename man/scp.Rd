\name{scp}
\alias{scp}
\title{Create a spatial community pattern}
\usage{
  scp(species, x = NULL, y = NULL, win = NULL,
    type = "ind-mapped", abundance = NULL, traits = NULL,
    habitat = NULL, phylo = NULL, forceUnique = FALSE)
}
\arguments{
  \item{species}{a vector of species names.}

  \item{abundance}{a vector of number of individuals for
  each species}

  \item{x,y}{spatial coordination of each individual if
  available}

  \item{win}{either a vector of length 4 or an owin object
  (\link{owin}) represents spatial region of the community}

  \item{traits}{a dataframe contains traits of every
  individual or species. e.g. DBH, status.}

  \item{habitat}{a list \link{im} objects about
  environments about the community. e.g. location, mean
  elevation and other topographic parameters of the
  quadrat.}

  \item{phylo}{a \link{phylo} object contains phylogeny of
  species in the community.}

  \item{type}{a character indicate types of the community.
  it can be "ind-mapped" (default), "pre-abs" and "sp-ab"
  which represents represents either a species
  presence-absence vector, species abundance vector or
  individual mapped dataframe of a community,
  respectively.}

  \item{forceUnique}{unique location of individuals if it
  set true.}
}
\value{
  a \code{scp} object contains following attributes:
  \tabular{ll}{ com: \tab a \link{ppp} object contains
  spatial information of the community \cr N: \tab total
  number of individual \cr S: \tab total number of species
  \cr sp: \tab a vector of total species name \cr ab: \tab
  species abundance distribution \cr win: \tab an
  \link{owin} object contains spatial region of the
  community \cr traits: \tab a dataframe contains character
  of each individuals, e.g. species name, DBH, status \cr
  habitat: \tab a list \link{im} objects about environments
  about the community. \cr phylo: \tab a \link{phylo}
  object contains phylogeny of species in the community.
  \cr }
}
\description{
  \code{scp} create a spatial community pattern object. it
  represents either a species presence-absence vector,
  species abundance vector or individual mapped dataframe
  of a community.
}
\details{
  if type="ind-mapped", the required parameters are
  \code{species}, \code{x}, \code{y} and \code{owin}.
  unique location check and removing duplicated indiviudals
  were carried out subsequently.
}
\examples{
#let's set a community has three species and 50 indiviudals with name A,B and C.
species=sample(c("A","B","C"),50,replace=TRUE)

#create a species presence-absence community
scp(species,type="pre-abs")

#create a species abundance community
scp(species,type="sp-ab")

#create a individual mapped community
x=runif(50,0,100)
y=runif(50,0,100)
win=owin(c(0,100),c(0,100))
scp(species,x=x,y=y,win=win,type="ind-mapped")

#an real scp object from BCI data
data(BCI)
BCI
}

