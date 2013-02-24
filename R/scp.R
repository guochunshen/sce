#' Create a spatial community pattern
#' 
#' \code{scp} create a spatial community pattern object. it represents either a species presence-absence vector,
#' species abundance vector or individual mapped dataframe of a community.
#' 
#' @param species a vector of species names.
#' @param abundance a vector of number of individuals for each species
#' @param x,y spatial coordination of each individual if available
#' @param win either a vector of length 4 or an owin object (\link{owin}) represents spatial region of the community
#' @param traits a dataframe contains traits of every individual or species. e.g. DBH, status.
#' @param habitat a list \link{im} objects about environments about the community. e.g. location, mean elevation and other topographic parameters of the quadrat.
#' @param phylo a \link{phylo} object contains phylogeny of species in the community.
#' @param type a character indicate types of the community. it can be "ind-mapped" (default), "pre-abs" and "sp-ab" which represents
#' represents either a species presence-absence vector, species abundance vector or individual mapped dataframe of a community, respectively.
#' @param forceUnique unique location of individuals if it set true.
#' 
#' @details
#' if type="ind-mapped", the required parameters are \code{species}, \code{x}, \code{y} and \code{owin}.
#' unique location check and removing duplicated indiviudals were carried out subsequently.
#' 
#' 
#' 
#'
#' @return a \code{scp} object contains following attributes:
#' \tabular{ll}{
#' com: \tab a \link{ppp} object contains spatial information of the community \cr
#' N: \tab  total number of individual \cr
#' S: \tab  total number of species \cr
#' sp: \tab  a vector of total species name \cr
#' ab: \tab  species abundance distribution \cr
#' win: \tab an \link{owin} object contains spatial region of the community \cr
#' traits: \tab a dataframe contains character of each individuals, e.g. species name, DBH, status \cr
#' habitat: \tab a list \link{im} objects about environments about the community. \cr
#' phylo: \tab  a \link{phylo} object contains phylogeny of species in the community. \cr
#' }
#' 
#' 
#' 
#' @examples
#' #let's set a community has three species and 50 indiviudals with name A,B and C.
#' species=sample(c("A","B","C"),50,replace=TRUE)
#' 
#' #create a species presence-absence community
#' scp(species,type="pre-abs")
#' 
#' #create a species abundance community
#' scp(species,type="sp-ab") 
#' 
#' #create a individual mapped community
#' x=runif(50,0,100)
#' y=runif(50,0,100)
#' win=owin(c(0,100),c(0,100))
#' scp(species,x=x,y=y,win=win,type="ind-mapped")
#' 
#' #an real scp object from BCI data
#' data(BCI)
#' BCI
#' 

scp <- function(species,x=NULL,y=NULL,win=NULL,type="ind-mapped",
                abundance=NULL,traits=NULL,habitat=NULL,phylo=NULL,forceUnique=FALSE){
  data=list()
  
  if("ind-mapped"==type){
    valid=check_ind_mapped_data(species=species, x=x, y=y, win=win, forceUnique=forceUnique)
    if(!valid) return
    com=data.frame(species=species,x=x,y=y)
    del_index=attr(valid,"del_index")
    if(forceUnique & !is.null(del_index)){
      com=com[-del_index,]
      traits=traits[-del_index,]
      species=species[-del_index]
    }
      
    data$com=ppp(x=com$x,y=com$y,window=win,check=FALSE)
    data$win=win
    data$N=dim(com)[1]
    data$ab=table(com$species)
    data$sp=unique(com$species)
    data$S=length(data$sp)
  }else if("pre-abs"==type){
    #TODO 
  }else if("sp-ab"==type){
    #TODO
  }else{
    stop("unsupported data types")
  }
  data$type=type
  traits$species=species
  data$traits=traits
  data$habitat=habitat
  
  class(data)="scp"
  return(data)
}


check_ind_mapped_data<-function(species,x=NULL,y=NULL,win=NULL,forceUnique=FALSE){
  valid=TRUE
  #check required paramters
  if(is.null(species) | is.null(x) | is.null(y) | is.null(win)){
    stop("\n All of the species, x, y and win parameters are needed.\n")
    valid=FALSE
    return
  }
  
  #check types of win
  if(!is(win,"owin")){
    stop("class of win is not owin")
  }
  
  #check NA values in species
  if(any(is.na(species)))
    warning("NA values are contained in species")
  
  #check finite of x and y
  if(any(!is.finite(x)) | any(!is.finite(y))){
    stop("\n infinite/NA values are not alowed in the coordination of individuals")
    valid=FALSE
    return
  }
  
  #check types of the arguments
  if(!is(x,"numeric") | !is(y,"numeric")){
    stop("Types of x and y should be numeric.\n")
    valid=FALSE
    return
  }
  
  #check length of the species, x, y
  if(length(unique(c(length(x),length(y),length(species))))!=1){
    stop("lengthes of species, x and y are not equal.")
  }
  
  #check uniqueness, maximum resolution of numeric value is six.
  x=round(x,6)
  y=round(y,6)
  xy=paste(x,y)
  xyfrq=table(xy)
  if(any(xyfrq>1)){
    if(!forceUnique){
      stop("locations of individual are not unique.\n")
    }else{
      del_index=match(names(xyfrq)[which(xyfrq>1)],xy)
      attr(valid,"del_index")=del_index
    }
  }
  
  #check the 
  return(valid)
}

#display main information of the scp object 
print.scp <- function(x,...){
  if("ind-mapped"==x$type){
    cat("A fully individual mapped community \n ")
    cat("with",x$N,"individuals belong to",x$S,"species ")
    cat("in an area of",area.owin(x$win))
  }
}



