#'
#' PCA transformation of the habitat variables
#'
#'@param habitat a list of im with the same dim
#'@param scale scale all variables under the unit variance.
#'


pcaHabitat=function(habitat,scale=TRUE){
  n=length(habitat)
  da=data.frame(as.vector(habitat[[1]]$v))
  for(i in 2:n){
    da=cbind(da,as.vector(habitat[[i]]$v))
  }
  pcaall=scores(rda(da,scale=scale),choices=1:n)$sites
  covr=dataframeToIms(pcaall,dim(habitat[[1]]),habitat[[1]]$xcol,habitat[[1]]$yrow)
  return(covr)
}