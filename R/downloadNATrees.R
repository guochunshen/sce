#'
#'Download theMaps of the ranges of tree species in North America compiled by Elbert Little, of the U.S. Department of Agriculture, Forest Service, and others (see references below) were digitized for use in USGS' vegetation-climate modeling studies.
#'
#'@param PATH the location where to save the data
#' 
#'@details
#' shapefiles at:
#'  http://esp.cr.usgs.gov/data/atlas/litle/
#'
#' zipped shapefile URL looks like:
#'  http://esp.cr.usgs.gov/data/atlas/little/abieamab.zip
#'  
#' Spatial_Reference_Information:
#' Horizontal_Coordinate_System_Definition:
#' Geographic:
#' Latitude_Resolution: 0.01
#' Longitude_Resolution: 0.01
#' Geographic_Coordinate_Units: Decimal degrees
#' Geodetic_Model:
#' Horizontal_Datum_Name: North American Datum of 1927
#' Ellipsoid_Name: Clarke 1866
#' Semi-major_Axis: 6378206.4
#' Denominator_of_Flattening_Ratio: 294.98
#' 
#' Spatial_Domain:
#' Bounding_Coordinates:
#' West_Bounding_Coordinate: -180.0
#' East_Bounding_Coordinate: -50.0
#' North_Bounding_Coordinate: 83.0
#' South_Bounding_Coordinate: 7.0
#' 
#' @author Jian Zhang "zjianATualberta.ca"
#'
#'@seealso
#'\code{\link{NAtreeRaster}}
#'

downloadNATree <- function(PATH="./"){
  REMOTE=TRUE
  INDEX_URL <- 'http://esp.cr.usgs.gov/data/atlas/little/'
  
  if(REMOTE){
    cat('checking index page...\n')
    sp.page <- scan(INDEX_URL, sep = '\n', what = '')
    start.data <- grep('Latin Name', sp.page) + 1
    
    sp.codes <- sp.page[start.data:length(sp.page)]
    sp.codes <- gsub('^[^"]*', '', sp.codes,
                     perl = TRUE, useBytes = TRUE)
    sp.codes <- gsub('.pdf.*$', '', sp.codes,
                     perl = TRUE, useBytes = TRUE)
    sp.codes <- gsub('\"', '', sp.codes,
                     perl = TRUE, useBytes = TRUE)
    sp.codes <- sp.codes[!sp.codes == '']
    sp.codes <- sp.codes[1:(length(sp.codes)-2)]
    #  sp.codes <- sp.codes[1:5] # testing
    
    if(length(sp.codes) != length(unique(sp.codes))){
      warning('there are non-unique 8-letter species codes')
    }
    cat(paste('there are', length(sp.codes),
              'species listed on the index page\n'))
    
    # Read in a hand-corrected species table....
    #  sp.table <- read.table(paste(PATH, '09-02-10LittleSpp.csv',
    #                         sep = ''), sep = '\t', header = TRUE)
    
    # there are 679 spp. as of 09-03-13
    
    # Download zipped shape files (with index page)
    oldwd <- getwd()
    setwd(PATH)
    if(!'zipped_shpfiles' %in% dir()){
      system('mkdir zipped_shpfiles')
    }
    cat('downloading index page\n')
    #  system(paste('/sw/bin/wget -N ', INDEX_URL, sep = ''))
    #  setwd(paste(PATH, 'zipped_shpfiles', sep = ''))
    for(i in 1:length(sp.codes)){
      print(i)
      filename <- paste(sp.codes[i], '.zip', sep = '')
      if(!filename %in% dir()){
        cat(paste('downloading zipped shapefile ',
                  sp.codes[i], '...\n', sep = ''))
        download.file(paste(INDEX_URL,filename,sep=""), destfile = filename) 				
        # try(system(paste('/sw/bin/wget ', INDEX_URL, filename, sep = '')))
      }
    }
    zippedfiles <- dir()
    setwd(oldwd)
  }else{
    zippedfiles <- dir(paste(PATH, 'zipped_shpfiles', sep = ''))
  }
  
  cat(paste('there are', length(zippedfiles),
            'zipped shape files in the right place\n'))
  
  
}