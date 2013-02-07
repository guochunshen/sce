#' Remove rare species from the community
#' 
#' @param com a scp object representing community
#' @param maxN the maximum number of individuals to be a rare species
#' 

removeRareSpecies<-function(com,maxN=10){
  rare_index=which(com$ab<=maxN)
  if(length(rare_index)==0)
    return(com)
  del_flag=rep(FALSE,com$N)
  for(i in rare_index){
    onesp=com$sp[i]
    del_flag=del_flag | (com$traits$species==onesp)
  }
  
  com$com=com$com[!del_flag]
  com$S=com$S-length(rare_index)
  com$N=com$N-sum(del_flag)
  com$traits=com$traits[!del_flag,]
  com$traits$species=as.factor(as.character(com$traits$species))
  com$ab=table(com$traits$species)
  com$sp=unique(com$traits$species)
  return(com)
}