#' Select species with specific abundances (>minN & <maxN) from the community
#' 
#' @usage
#' selSpeciesByAbund(com,minN,maxN)
#' removeRareSpecies(com,minN=10)
#' 
#' @param com a scp object representing community
#' @param maxN the maximum number of individuals for a species, not included
#' @param minN the minimum number of individuals for a species, not included
#' 
#' 
#' @examples
#' 
#' data(testData)
#' 
#' com=selSpeciesByAbund(testData,minN=10,maxN=50)
#' 
#' range(com$ab) #11 49
#' 
#' com2=removeRareSpecies(testData,minN=5)
#' 
#' min(com2$ab)>5
#' 
#' 

selSpeciesByAbund<-function(com,minN=10,maxN=50){
  del_index=which(com$ab<=minN | com$ab>=maxN)
  if(length(del_index)==0){
    return(com)
  }
  del_flag=rep(FALSE,com$N)
  for(i in del_index){
    onesp=com$sp[i]
    del_flag=del_flag | (com$traits$species==onesp)
  }
  com=subset(com,!del_flag)
  return(com)
}


removeRareSpecies<-function(com,minN=10){
  return(selSpeciesByAbund(com,minN=minN,maxN=Inf))
}