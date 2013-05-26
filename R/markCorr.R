#'
#' The c++ implemented mark correlation function with mark random labeling null model
#'
#'@param com a scp object
#'@param markName the name of the mark contained in \code{com} to calculate the mark pattern
#'@param r a vector of neighborhood distances to quantifying the mark correlation function
#'@param testFunName name of the test function, see details of the current implemented test functions.
#'@param nsim number of shuffling the values of mark
#'@param h the torrent bin.
#'@param exclude_conspecific a logical flag to whether include the conspecific indiviual or not
#'@param normalize a logical flag to normalize the mark correlation function
#'@param alpha significant level used to calculate pointwise confidence interval
#'
#'@details
#'it calculate the mark correlation function for individual based numerical marks, e.g. DBH.
#'current implemented test functions are
#'
#'testFunName: abdif (absolute difference): abs(m1-m2)
#'testFunName: diff (difference from focal to neighbor): m1-m2
#'testFunName: sum (sum of the marks of two points): m1+m2
#'testFunName: reldif (relative absolute difference): abs(m1-m2)/(m1+m2)
#'
#'
#'@examples
#'data(testData)
#'re1=markCorr(com=testData,markName="dbh",r=1:10,testFunName="tf1",nsim=10,h=1,exclude_conspecific=TRUE)
#'
#'

markCorr<-function(com,markName,r,testFunName,nsim=10,h=0.5,exclude_conspecific=FALSE,
                   normalize=TRUE,alpha=0.05){
  mi=which(names(com$traits)==markName)
  if(length(mi)==0)
    stop("no mark is found equals to the given mark name")
  mark=com$traits[,mi]
  if(inherits(mark,"logical") | inherits(mark, "character"))
    stop("The mark is not numeric, please use the mark connection function instead")
  
  #remove missing marks
  deli=is.na(mark)
  if(sum(deli)>0){
    mark=mark[!deli]
    com=subset(com,!deli)
  }
  #number of total individuals
  N=length(mark)
  randomMarki=replicate(nsim,sample(1:N))
  
  sp=as.numeric(as.factor(com$traits$species))
  
  tfnames=c("abdif","diff","sum","reldif")
  
  tftype=which(tfnames==testFunName)-1
  
  #void markcorr(int* N, double* x, double* y, double* mark, double* r, int* nr, double* h, int* tftype,
  #int* nsim, int* marki, int* exclude_conspecific,double* mkcvalues)
  nr=length(r)
  marki=c(1:N,as.vector(randomMarki))
  mkcvalues=rep(0,times=nr*(nsim+1))
  const=0;
  
  re=.C("markcorr",as.integer(N),as.double(com$com$x),as.double(com$com$y),as.double(mark),as.double(r),
     as.integer(nr),as.double(h),as.integer(tftype),as.integer(nsim),as.integer(marki),
     as.integer(exclude_conspecific),as.integer(sp),as.double(mkcvalues),as.integer(normalize),
        as.double(const))[c(13,15)]
  
  const=re[[2]]
  re=re[[1]]
  dim(re)=c(nr,nsim+1)
  mkcobs=re[,1]
  if(nsim!=0){
    mkcsim=re[,-1]
    conf_low_index=floor(nsim*alpha/2)
    if(conf_low_index==0)
      conf_low_index=1
    conf=apply(mkcsim,1,function(x) sort(x)[c(conf_low_index,nsim-conf_low_index)])
    dim(conf)=c(2,nr)
    return(list(obs=mkcobs,low=conf[1,],upper=conf[2,],const=const))
  }else{
    return(list(obs=mkcobs,const=const))
  }
  
}

