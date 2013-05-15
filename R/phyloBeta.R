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
#'@examples
#'library(testthat)
#'
#'com=rCom(1000,10,win=owin(c(0,100),c(0,100)),ab="physignal",phy=list(br=rexp,phylosignal=1000))
#'com=quadratize(com,20,20)
#'phyd=cophenetic(com$phylo)
#'re1=phyloBeta(com,phyd=phyd,Fun=comdist_C,nsim=10,abundance.weighted=FALSE)
#'re2=phyloBeta(com,phyd=phyd,Fun=comdist_C,nsim=10,abundance.weighted=TRUE)
#'
#'expect_true(any(re1$real!=re2$real))
#'
#'re1=phyloBeta(com,phyd=phyd,Fun=comdistnt_C,nsim=10,abundance.weighted=FALSE)
#'re2=phyloBeta(com,phyd=phyd,Fun=comdistnt_C,nsim=10,abundance.weighted=TRUE)
#'
#'expect_true(any(re1$real!=re2$real))
#'
#'expect_true(max(re1$r)<50)
#'
#'

phyloBeta<-function(com,phyd,Fun,nsim,rmax=50,alpha=0.05,...){
  #get the community matrix
  comtable=table(com$traits$ploti,com$traits$species)
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
  
  phylo_beta_nulls=unlist(lapply(phy_nulls,function(x) {
    tapply(as.numeric(Fun(comm=comtable,dis=x,cal_pairs=cal_pairs, ...)),spaced,mean,na.rm=TRUE)
    
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
