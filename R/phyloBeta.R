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

phyloBeta<-function(com,phyd,Fun,nsim,rmax=50,alpha=0.05,...){
  #get the community matrix
  comtable=table(com$traits$ploti,com$traits$species)
  #get the community spatial distance matrix
  quadratindex=as.numeric(rownames(comtable))
  spaced=as.matrix(dist(attr(com,"quadratxy")[quadratindex,]))
  spaced=spaced[lower.tri(spaced)]
  
  #calcualte the observed comdist
  phylo_beta_obs=tapply(as.numeric(Fun(comm=comtable,dis=phyd, ...)),spaced,mean,na.rm=TRUE)
  
  #using species shuffling null model to generate a phylogeny tree under the null model
  phy_nulls=spShuffle(phyd,nsim)
  
  #calculate the confidence envelope
  phylo_beta_nulls=unlist(lapply(phy_nulls,function(x) {
    tapply(as.numeric(Fun(comm=comtable,dis=x, ...)),spaced,mean,na.rm=TRUE)
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
  selr=r<=rmax
  re=list(r=r[selr],real=as.numeric(phylo_beta_obs)[selr],
          lower=phylo_beta_conf[1,][selr],upper=phylo_beta_conf[2,][selr],pvalues=pvalues[selr])
  return(re)
}