\name{selSpeciesByAbund}
\alias{selSpeciesByAbund}
\title{Select species with specific abundances (>minN & <maxN) from the community}
\usage{
  selSpeciesByAbund(com,minN,maxN)
    removeRareSpecies(com,minN=10)
}
\arguments{
  \item{com}{a scp object representing community}

  \item{maxN}{the maximum number of individuals for a
  species, not included}

  \item{minN}{the minimum number of individuals for a
  species, not included}
}
\description{
  Select species with specific abundances (>minN & <maxN)
  from the community
}
\examples{
data(testData)

com=selSpeciesByAbund(testData,minN=10,maxN=50)

range(com$ab) #11 49

com2=removeRareSpecies(testData,minN=5)

min(com2$ab)>5
}

