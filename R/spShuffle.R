#'
#'shuffling tip label among tips of phylogeny
#'
#'@param phyd a phylogenetic distance matrix with col and row species names
#'@param n number of shuffling
#'
#'
#'@return 
#'a list of length n about the shuffled phylogenetic distance matrix
#'

spShuffle<-function(phyd,n){
  spnames=colnames(phyd)
  re=list()
  for(i in 1:n){
    sprandom=sample(spnames)
    rownames(phyd)=colnames(phyd)=sprandom
    re[[i]]=phyd
  }
  return(re)
}