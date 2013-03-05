\name{NAtreeRaster}
\alias{NAtreeRaster}
\title{transfer the zipped shapefile of species distribution into a species occurance data.frame}
\usage{
  NAtreeRaster(PATH, gridsize = 0.5)
}
\arguments{
  \item{PATH}{the location of zip file of species
  distribution}

  \item{gridsize}{the size of grid to divide map of north
  american (range: -180, -50, 7, 83).}
}
\value{
  a dataframe, the first and second columuns are the
  spatial location of the grid, other columuns represent
  distribution of a species. NA means not observed that
  species in that grid.
}
\description{
  transfer the zipped shapefile of species distribution
  into a species occurance data.frame
}
\author{
  Jian Zhang "zjianATualberta.ca"
}
\seealso{
  \code{\link{downloadNATree}}
}

