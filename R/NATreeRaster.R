#'
#'transfer the zipped shapefile of species distribution into a species occurance data.frame
#'
#'@param PATH the location of zip file of species distribution 
#'@param gridsize the size of grid to divide map of north american (range: -180, -50, 7, 83). 
#'
#'@return
#'a dataframe, the first and second columuns are the spatial location of the grid, other columuns represent distribution 
#'of a species. NA means not observed that species in that grid.
#'
#' @author Jian Zhang "zjianATualberta.ca"
#' 
#'@seealso
#'\code{\link{downloadNATree}}

NAtreeRaster<-function(PATH,gridsize=0.5){
  prePath=getwd()
  setwd(PATH)
  zipfiles <- dir() # show all zip files
  zipfiles
  
  # raster
  library(raster)
  r <- raster(extent(c(-180, -50, 7, 83)))
  # GridSize (unit: degree)
  res(r)=gridsize
  xy <- xyFromCell(r,1:ncell(r)) # LongLat
  
  spOccur <- matrix(NA,dim(xy)[1],length(zipfiles)+2)
  spOccur[,1] <- xy[,1] # long
  spOccur[,2] <- xy[,2] # lat
  
  for (i in 1:length(zipfiles)){
 
    # unzip 
    unzip(zipfiles[i])
    
    filename <- strsplit(zipfiles[i],".zip")[[1]]
    
    # Read shapefile into R
    library(rgdal)
    shpfile <- try(readOGR(getwd(),filename))
    if(inherits(shpfile,"try-error")){
      warning(paste("The",i,"species occured some error!"))
      next()
    }
      
    extent(shpfile)
    tst <- try(rasterize(shpfile, field=names(shpfile@data)[3], r))
    if(inherits(tst,"try-error")){
      warning(paste("The",i,"species occured some error!"))
      next()
    }
      
    #tst <- rasterize(shpfile, field=paste(toupper(filename),"_",sep=""), r)
    # tst <- rasterize(shpfile, r, fun="count")
    tst[tst>0] <- 1
    
    print(i)
    print(range(tst[],na.rm=T))
    
    # mapping
    #plot(tst,col="red",main=filename)
    #library(maps)
    #map("world",xlim=c(-180,-50),ylim=c(7,83),add=TRUE)
    
    # convert to data.frame
    spOccur[,2+i] <- tst[]
  }
  
  rm(i,xy,tst,r,shpfile)
  
  # as.data.frame
  spOccur <- as.data.frame(spOccur)
  colnames(spOccur) <- c("long","lat",do.call("rbind",strsplit(zipfiles,".zip"))[,1])
  setwd(prePath)
  return(spOccur)
}