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
#' subset(testData,1:100)
#' 
#' #get a perticular species
#' sel=testData$traits$species==testData$sp[2]
#' subset(testData,sel)
#' 



"[.scp" <-
  function(x, i,...) {
    
  #todo some times it works, some times are not, figure it out   
  x$com=x$com[i]
  x$N=x$com$n
  x$traits=x$traits[i,]
  x$traits$species=as.factor(as.character(x$traits$species))
  x$ab=table(x$traits$species)
  x$sp=unique(x$traits$species)
  x$S=length(x$sp)
  return(x)
}

eval("[.scp")


subset.scp <- function(x, i){
  
  x$com=x$com[i]
  x$N=x$com$n
  x$traits=x$traits[i,]
  x$traits$species=as.factor(as.character(x$traits$species))
  x$ab=table(x$traits$species)
  x$sp=unique(x$traits$species)
  x$S=length(x$sp)
  
  return(x)
}

