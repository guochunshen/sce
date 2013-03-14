#'
#' Remvoe species from the community by a name vector
#' 
#'@param scp an scp object
#'@param species_names a vector of species name that will be removed from the community
#'
#'@return
#'a scp object without the given species
#'
removeSpecies<-function(scp,species_names){
  ns=length(species_names)
  if(ns==0)
    return(scp)
  delflag=rep(FALSE,scp$N)
  for(i in 1:ns){
    delflag=delflag | scp$traits$species == species_names[i]
  }
  return(subset(scp,!delflag))
}