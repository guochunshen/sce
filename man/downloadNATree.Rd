\name{downloadNATree}
\alias{downloadNATree}
\title{Download theMaps of the ranges of tree species in North America compiled by Elbert Little, of the U.S. Department of Agriculture, Forest Service, and others (see references below) were digitized for use in USGS' vegetation-climate modeling studies.}
\usage{
  downloadNATree(PATH = "./")
}
\arguments{
  \item{PATH}{the location where to save the data}
}
\description{
  Download theMaps of the ranges of tree species in North
  America compiled by Elbert Little, of the U.S. Department
  of Agriculture, Forest Service, and others (see
  references below) were digitized for use in USGS'
  vegetation-climate modeling studies.
}
\details{
  shapefiles at: http://esp.cr.usgs.gov/data/atlas/litle/

  zipped shapefile URL looks like:
  http://esp.cr.usgs.gov/data/atlas/little/abieamab.zip

  Spatial_Reference_Information:
  Horizontal_Coordinate_System_Definition: Geographic:
  Latitude_Resolution: 0.01 Longitude_Resolution: 0.01
  Geographic_Coordinate_Units: Decimal degrees
  Geodetic_Model: Horizontal_Datum_Name: North American
  Datum of 1927 Ellipsoid_Name: Clarke 1866
  Semi-major_Axis: 6378206.4
  Denominator_of_Flattening_Ratio: 294.98

  Spatial_Domain: Bounding_Coordinates:
  West_Bounding_Coordinate: -180.0
  East_Bounding_Coordinate: -50.0
  North_Bounding_Coordinate: 83.0
  South_Bounding_Coordinate: 7.0
}
\author{
  Jian Zhang "zjianATualberta.ca"
}
\seealso{
  \code{\link{NAtreeRaster}}
}

