#' Extract a subject set of a community
#' 
#' @param x an scp object
#' @param i in numeric or logical index of needed data
#' 
#' 
#' 
#' @examples
#' data(testData)
#' 
#' #get the first one hundred individuals
#' subsetScp(testData,1:100)
#' 
#' #get a perticular species
#' sel=testData$traits$species==testData$sp[2]
#' subsetScp(testData,sel)
#' 

subsetScp <- function(x, i){
  x$com=x$com[i]
  x$N=x$com$n
  x$traits=x$traits[i,]
  x$traits$species=as.factor(as.character(x$traits$species))
  x$ab=table(x$traits$species)
  x$sp=unique(x$traits$species)
  x$S=length(x$sp)
  return(x)
}
 