#'
#'calculate the the quadrat based phylogenetic beta diversity with significant test by simulations
#'
#'@param com a quadratized scp object trait for each individual and phylogeney
#'@param phyd phylogenetic distance matrix
#'@param nsim number of shuffling tip from the given phylongeny
#'@param rmax maximum spatial distance between two communities
#'@param ... other parameters past to \code{\link{comdist_C}} function.
#'
#'
#'
#'
#'

phyloBeta<-function(com,phyd,Fun,nsim,rmax=50,alpha=0.05,...){
  #get the community matrix
  comtable=table(com$traits$ploti,com$traits$species)
  #if there is too much quadrat that needed to be calculate the pairwise dissimilarity, using a different stratage
  if(dim(comtable)[1]>15000){
    return(phyloBeta2(com,phyd,Fun,nsim,rmax,alpha))
  }
  #get the community spatial distance matrix
  quadratindex=as.numeric(rownames(comtable))
  spaced=as.matrix(dist(attr(com,"quadratxy")[quadratindex,]))
  #define the community pairs within the maximum distance away
  cal_pairs=spaced<=rmax
  spaced[!cal_pairs]=NA
  spaced=round(spaced[lower.tri(spaced)],2)
  
  
  #calcualte the observed comdist
  phylo_beta_obs=tapply(as.numeric(Fun(comm=comtable,dis=phyd, cal_pairs=cal_pairs, ...)),spaced,mean,na.rm=TRUE)
  
  #using species shuffling null model to generate a phylogeny tree under the null model
  phy_nulls=spShuffle(phyd,nsim)
  
  #calculate the confidence envelope
  
  phylo_beta_nulls=unlist(lapply(1:nsim,function(x) {
    print(x)
    tapply(as.numeric(Fun(comm=comtable,dis=phy_nulls[[x]],cal_pairs=cal_pairs, ...)),spaced,mean,na.rm=TRUE)
    
  }))
  dim(phylo_beta_nulls)=c(length(phylo_beta_obs),nsim)
  lower=round(nsim*alpha/2)
  if(lower==0)
    lower=1
  upper=round(nsim*(1-alpha/2))
  phylo_beta_conf=apply(phylo_beta_nulls,1,function(x) {x=sort(x); x=x[!is.na(x)];x[c(lower,upper)]})
  dim(phylo_beta_conf)=c(2,length(phylo_beta_obs))
  
  #calculate the pvalue
  position=phylo_beta_nulls>as.numeric(phylo_beta_obs)
  pvalues=apply(position,1,sum)/nsim
  pvalues[pvalues>0.5]=1-pvalues[pvalues>0.5]
  
  r=as.numeric(names(phylo_beta_obs))

  re=list(r=r,real=as.numeric(phylo_beta_obs),
          lower=phylo_beta_conf[1,],upper=phylo_beta_conf[2,],pvalues=pvalues)
  return(re)
}

#a different interface
phyloBeta2<-function(com,phyd,type="dpw",nsim,rmax=50,alpha=0.05,abundance.weighted=TRUE,exclude_conspecifics=TRUE,...){
  
  comtable=table(com$traits$ploti,com$traits$species)
  qxy=as.matrix(attr(com,"quadratxy"))

  phylo_beta_obs=phyloBeta_one(comtable, phyd, type, abundance.weighted, qxy, rmax, exclude_conspecifics)
  
  phy_nulls=spShuffle(phyd,nsim)
  phylo_beta_nulls=unlist(lapply(phy_nulls,function(x) {
   
    re=phyloBeta_one(comtable, x, type, abundance.weighted, qxy, rmax, exclude_conspecifics)
    return(re)
  }))
  
  dim(phylo_beta_nulls)=c(length(phylo_beta_obs),nsim)
  lower=round(nsim*alpha/2)
  if(lower==0)
    lower=1
  upper=round(nsim*(1-alpha/2))
  phylo_beta_conf=apply(phylo_beta_nulls,1,function(x) {x=sort(x); x=x[!is.na(x)];x[c(lower,upper)]})
  dim(phylo_beta_conf)=c(2,length(phylo_beta_obs))
  
  #calculate the pvalue
  position=phylo_beta_nulls>as.numeric(phylo_beta_obs)
  pvalues=apply(position,1,sum)/nsim
  pvalues[pvalues>0.5]=1-pvalues[pvalues>0.5]
  
  r=as.numeric(names(phylo_beta_obs))
  
  re=list(r=r,real=as.numeric(phylo_beta_obs),
          lower=phylo_beta_conf[1,],upper=phylo_beta_conf[2,],pvalues=pvalues)
  return(re)
}

phyloBeta_one <- function (comtable, phyd, type, abundance.weighted, qxy, rmax, exclude_conspecifics) {
  dat <- match.comm.dist(as.matrix(comtable), phyd)
  x <- dat$comm
  dis <- as.matrix(dat$dist)
  
  
  if(type=="dpw"){
    indextype=0
    if (!abundance.weighted) {
      x <- decostand(x, method = "pa")
    }
    x <- decostand(x, method = "total", MARGIN = 1)
  } else if(type=="dnn"){
    indextype=1
    x <- decostand(x, method = "total", MARGIN = 1)
    if(abundance.weighted){
      x=t(apply(x,1,function(x) x/sum(x)/2))
    }else{
      x[x!=0]=1
      x=x/apply(x,1,sum)/2
    }
  }
  
  N <- dim(x)[1]
  S <- dim(x)[2]  
  #print(list(N, dis, x,qxy,rmax,exclude_conspecifics, indextype))
  #observe
  phyobs=phyloBeta_C(N, dis, x,qxy,rmax,exclude_conspecifics, indextype)
  phylo_beta_obs=tapply(phyobs$phyd,round(phyobs$spaced,2),mean,na.rm=TRUE)
  return(phylo_beta_obs)
}