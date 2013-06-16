#' divide the whole community into small and large DBH groups
#' 
#' @param com a scp object
#' 
#' @return
#' a logical vector of small individual which is true
#' 
#' 

smallLargeDBH=function(com){
  smallij=rep(FALSE,com$N)
  for(sp in com$sp){
    spi=which(com$traits$species==sp)
    dbh=com$traits$dbh[spi]
    small=(dbh<=quantile(dbh,na.rm=TRUE)[3])
    smallij[spi]=small
  }
  return(smallij)
}