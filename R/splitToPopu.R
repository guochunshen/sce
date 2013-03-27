#'
#'Split the whole community into a list of populations
#'
#'@param com a scp object
#'

splitToPopu<-function(com){
  splist=com$sp
  re=list()
  for(i in 1:com$S){
    re[[i]]=subset(com,com$traits$species==splist[i])
  }
  names(re)=splist
  return(re)
}