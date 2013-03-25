#'
#'create a community driven by single process like dispersal limitation, habitat association or competition
#'
#'@param com_control a list object contains expected properties about the simulated community.
#'
#'@details
#'
#'The com_control has follow field about the community:
#'
#'N: total number individuals in the community
#'
#'S: total number species in the community
#'
#'win: an \code{\link{owin}} object represent the spatial range of the community
#'
#'ab: a character discrebing the species abundance distribution of the target community. currently avaible options are "uniform", "logseries" and "lognormal"
#'
#'intra: a character representing intraspecific relationship. currently avaible options are "Poisson", "dispersal"
#'
#'inter: a string representing interspecific relationship. currently avaible options are "independent" and "dependent"
#'
#'covr: a logical flag to add species habitat association into the community. currently, habitat is simpled modeled by the sine function.
#'      it always be FALSE in community assembled by competition.
#'
#'covrscale: the spatial scale of the covariance.
#'
#'beta: a real number defines the competition intensity between species. the large the number is, the intense the competition is.
#'
#'r: a real distance value. competition only happend within this distance.
#'
#'R2: define the strength of correlation between species character and their phylogenetic relationship
#'
#'
#'@return
#'a list contains: a marked point pattern and a ecological difference matrix
#'
#'@examples
#'
#'  com_control=list(S=10,N=1000,ab="lognormal",win=owin(c(0,100),c(0,100)),
#' intra="Poisson",inter="independent",R2=0.9,covr=FALSE,beta=2,r=5,covrscale=16)
#' 
#' #Poisson distribution
#' rCommunity(com_control)
#'
#' #dispersal 
#' com_control$intra="dispersal"
#' rCommunity(com_control)
#' 
#' #habitat
#' com_control$intra="Poisson"
#' com_control$covr=TRUE
#' rCommunity(com_control)
#' 
#' #compeition, relatively time consuming
#' com_control$covr=FALSE
#' com_control$inter="dependent"
#' rCommunity(com_control)
#'


rCommunity=function(com_control){
  if(com_control$inter=="independent"){
    com=random_community(com_control)
  }else{
    com=try(complex_community(com_control))
    while(class(com)=="try-error"){
      com=try(complex_community(com_control))
    }
  }
  return(com)
}



random_community=function(com_control){
  N=com_control$N
  S=com_control$S
  win=com_control$win
  area=com_control$win$xrange[2]*com_control$win$yrange[2]
  #total intensity
  if(com_control$ab=="uniform")
    N_intensity=rep(N/S,S)/area
  else if(com_control$ab=="logseries"){
    N_intensity=rlogseries(N,S)/area
  }else if(com_control$ab=="lognormal"){
    N_intensity=rlognormal(N,S)/area
  }
  #species name
  spnames=1:com_control$S
  
  #random independent multiple point pattern
  if(com_control$intra=="Poisson"){
    if(com_control$covr){
      niche=matrix(c(runif(S,-1,1)),nrow=S,ncol=1)
      trend=list()
      covr=sin(seq(0,win$xrange[2])/com_control$covrscale)
      cl=length(seq(0,win$xrange[2]))
      covr=matrix(rep(covr,each=length(covr)),ncol=cl,nrow=cl)
      for(i in 1:S){
        tri=exp(abs(covr-niche[i]))
        
        trend[[i]]=im(tri,xrange=win$xrange,yrange=win$yrange)
      }
      data.ppp=rmpoint(com_control$N,trend,types=spnames,ptypes=N_intensity/sum(N_intensity,na.rm=T),win=com_control$win)
      attr(data.ppp,"covr")=covr
    }else{
      data.ppp=rmpoint(com_control$N,N_intensity,types=spnames,win=com_control$win)	
    }
    
    
  }else if(com_control$intra=="dispersal"){
    
    for(i in 1:S){
      N_exp=round(N_intensity[i]*area)
      kappa=N_intensity[i]/8
      X <- rThomas(kappa, 5, 8, win=com_control$win)
      if(X$n < N_exp & X$n > N_exp*1.2){
        X <- rThomas(kappa, 5, 8, win=com_control$win)
      }
      if(X$n>N_exp){
        delxy=sample(1:X$n,X$n-N_exp)
        X=X[-delxy]
      }
      
      if(i==1){
        x=X$x
        y=X$y
        sp=rep(i,X$n)
      }else{
        x=c(x,X$x)
        y=c(y,X$y)
        sp=c(sp,rep(i,X$n))
      }
    }
    data.ppp=ppp(x=x,y=y,marks=as.character(sp),window=com_control$win)
  }
  
  data.ppp$marks=data.frame(sp=marks(data.ppp))
  if(!com_control$covr){
    result=list(data.ppp=data.ppp,ecold=as.matrix(dist(1:com_control$S)))
  }else{
    result=list(data.ppp=data.ppp,ecold=as.matrix(dist(niche)))
  }
  return(result)
}


complex_community=function(com_control){
  S=com_control$S
  N=com_control$N
  win=com_control$win
  area=win$xrange[2]*win$yrange[2]
  if(com_control$ab=="uniform")
    ab=rep(N/S,S)
  else if(com_control$ab=="logseries"){
    ab=rlogseries(N,S)
  }else if(com_control$ab=="lognormal"){
    ab=rlognormal(N,S)
  }
  N_intensity=ab/win$xrange[2]/win$yrange[2]
  niche=matrix(c(runif(S,0,1)),nrow=S,ncol=1)
  #intensity of each type
  beta=rep(com_control$beta,S)
  #strength of interaction
  gamma=as.matrix(dist(niche))
  gamma=round(gamma/max(gamma),2)
  gamma=sqrt(gamma+0.5)
  gamma=gamma/max(gamma)
  #diag(gamma)=0.00
  #gamma=gamma[S:1,]
  rownames(gamma)=colnames(gamma)=1:S
  #radio of interaction
  r=matrix(com_control$r,nrow=S,ncol=S)
  mod08 <- list(cif="straussm",par=list(beta=beta,gamma=gamma,radii=r),
                w=win)
  if(FALSE){
    trend=list()
    covr=sin(seq(0,win$xrange[2])/10)
    cl=length(seq(0,win$xrange[2]))
    covr=matrix(rep(covr,each=length(covr)),ncol=cl,nrow=cl)
    for(i in 1:S){
      tri=exp(abs(covr-niche[i]))/area
      mu=log(N_intensity[i]/sum(exp(tri))*area)
      trend[[i]]=im(exp(tri+mu)*20,xrange=win$xrange,yrange=win$yrange)
    }
    mod08=c(mod08,trend=list(trend))
  }
  
  nr   <- 5e5
  nv  <- 1e5
  X <- rmh(model=mod08,start=list(n.start=N),
           control=list(ptypes=ab/N,nrep=nr,nverb=nv))
  X$marks=data.frame(sp=marks(X))
  return(list(data.ppp=X,ecold=gamma))	
}


nichesurface=function(niche,win){
  return(density(da))
}
