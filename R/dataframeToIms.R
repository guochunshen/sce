#'
#' convert a data.frame into a list of im objects
#' 
#' @param data a data.frame object contains observed variables
#' @param size size of the observed variable matrix
#' @param xcol ector of x coordinates for the pixel gid
#' @param yrow ector of y coordinates for the pixel gid
#' 
#' @return
#' a list of im objects
#'
#'

dataframeToIms<-function(data, size, xcol, yrow){
  covr = list()
  for (i in 1:dim(data)[2]) 
    covr[[i]] = im(matrix(data[,i], size[1], size[2]), xcol = xcol, yrow = yrow)
  names(covr)=colnames(data)
  return(covr)
}