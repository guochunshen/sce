#'
#' The phylogenetic mark correlation function with species shuffling null model
#'
#'@param com a scp object 
#'@param phyd a phylogenetic distance matrix with species colnames and colnames
#'@param focalij a logical vector indicate which individuals are focal individuals.
#'      if it is null, all individuals are the focal points.
#'@param rmax the maximum distance to calculate the phylogenetic mark correlatin function
#'@param step the step of distance from 0 to quantify spaital phylogenetic structure
#'@param nsim number of shuffling tip labels in the species shuffling null model
#'@param alpha sinificant level
#'@param scale a logical to indicate wether the phylogenetic mark correlation function should be scaled.
#'
#'@details
#'explained detailed further
#'
#'@examples
#'data(testData)
#'
#'phyMarkCorr(testData,testData$phylo,rmax=20)
#'


phyMarkCorr=function(com,phyd,focalij=NULL,rmax,step=1,nsim=19,alpha=0.05,scale=TRUE){
   
   if(is.null(focalij)){
     focal=1:(com$N)-1
     nfocal=com$N
   }else{
     focal=which(focalij)-1
     nfocal=length(focal)
   }
  sp=com$traits$species
  tips=rownames(phyd)
  sp=match(sp,tips)-1
  tips=1:length(tips)-1
  ntotal=com$N
  nsp=length(tips)
  
  nsteps = ceiling((rmax)/(step));
  pk=rep(0,nsteps)
  upper=rep(0,nsteps)
  lower=rep(0,nsteps)
  mean=rep(0,nsteps)
  pvalues=rep(0,nsteps)
  
  scale=as.numeric(scale)
  
  re=.C("phylocorr", as.integer(focal), as.integer(nfocal), as.double(com$com$x), as.double(com$com$y),
        as.integer(ntotal), as.integer(sp), as.integer(nsp), as.double(rmax), 
        as.double(as.numeric(phyd)), as.double(step), as.double(pk), as.integer(nsim),
        as.double(mean), as.double(upper), as.double(lower), as.double(pvalues),
        as.double(alpha), as.integer(scale));
  r=seq(step,rmax,step)
  result=list("real"=re[[11]],"upper"=re[[14]],"lower"=re[[15]],"r"=r,"pvalues"=re[[16]])
  class(result)="pmc"
  return(result)
}


plot.pmc<-function(pmc,...){
  plot(x=range(pmc$r),y=range(unlist(pmc[1:3])),type="n",...)
  polygon(x=c(pmc$r,rev(pmc$r)),y=c(pmc$lower,rev(pmc$upper)),col="grey",border="grey")
  lines(x=pmc$r,y=pmc$real,col=2)
}